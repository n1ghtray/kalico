import serial
import logging
import json
import time

from serial import SerialException

class Inkjet:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.reactor = self.printer.get_reactor()

        self.printer.register_event_handler("klippy:ready", self._handle_connect)
        self.gcode = self.printer.lookup_object("gcode")
        self.toolhead = None
        self.current_layer = None

        # Offsets and speeds
        self.inking_axis = config.get("inking_axis", "x")
        self.x_offset = config.getfloat("x_offset", 0)
        self.y_offset = config.getfloat("y_offset", 0)
        self.z_offset = config.getfloat("z_offset", 1)
        self.z_lift = config.getfloat("z_lift", 0.2)
        self.z_lift_speed = config.getfloat("z_lift_speed", 20)
        self.inking_speed = config.getfloat("inking_speed", 50)
        self.inking_accel = config.getfloat("inking_accel", 1000)
        self.inking_accel_margin = config.getfloat("inking_accel_margin", 10)
        self.inking_pass_height = config.getfloat("inking_pass_height", 5.4)
        self.inking_accel_delay = config.getfloat("inking_accel_delay", 10000)
        self.inking_columns_delay = config.getfloat("inking_columns_delay", 64)
        self.inking_start_layer = config.getint("inking_start_layer", 0)
        self.positioning_speed = config.getfloat("positioning_speed", 100)
        self.positioning_accel = config.getfloat("positioning_accel", 1000)
        self.purge_x_position = config.getfloat("purge_x_position")
        self.purge_y_position = config.getfloat("purge_y_position")
        self.purge_z_lift = config.getfloat("purge_z_lift")
        self.purge_cycles = config.getint("purge_cycles", 500)
        self.purge_delay = config.getint("purge_delay", 5)

        # State and data
        self.is_inking = False
        self.position_before_inking = None
        self.conn = None
        self.serial = config.get("serial")
        self.trigger_pin = config.get("trigger_pin", "test")
        self.pictures_base_path = config.get("pictures_base_path", "home/pi/printer_data/gcode_ink/")
        self.pictures_directory = None
        self.current_layer_data = None
        self.current_layer_passes = None
        self.current_pass = None
        self.print_mesh_min = None
        self.print_mesh_max = None

        self.x_starting_point = 0
        self.y_starting_point = 0
        self.pass_size = 0

        # Register G-code commands
        self.gcode.register_command("INKJET_CONNECT_BOARD", self.connect_board)
        self.gcode.register_command("INKJET_CONFIGURE_BOARD", self.configure_board)
        self.gcode.register_command("INKJET_DISCONNECT_BOARD", self.disconnect_board)
        self.gcode.register_command("INKJET_GET_BOUNDS", self.get_bounds)
        self.gcode.register_command("INKJET_PURGE", self.move_to_purge)
        self.gcode.register_command("INKJET_SET_DIRECTORY", self.set_directory)
        self.gcode.register_command("INKJET_SET_LAYER", self.set_layer)
        self.gcode.register_command("INKJET_PRINT", self.print_layer)
        self.gcode.register_command("FILL_DEBUG_INFO", self.fill_debug_info)
        self.gcode.register_command("INKJET_FILL_TEST", self.fill_test)

    def _handle_connect(self):
        self.toolhead = self.printer.lookup_object("toolhead")
        return

    def fill_test(self, gcmd):
        self.write_layer_pass_to_board(gcmd, self.current_layer_data[0])
        return


    def fill_debug_info(self, gcmd):
        self.connect_board(gcmd)
        self.configure_board(gcmd)
        self.pictures_directory = "/home/pi/printer_data/gcode_ink/TEST"
        self.current_layer = gcmd.get_int("LAYER", 300)
        self.load_layer_data(gcmd)
        xmin = gcmd.get_int("XMIN", 10)
        xmax = gcmd.get_int("XMAX", 97)        
        ymin = gcmd.get_int("YMIN", 10)
        ymax = gcmd.get_int("YMAX", 65)
        self.print_mesh_min = [xmin,ymin]
        self.print_mesh_max = [xmax,ymax]
        self.x_starting_point = self.print_mesh_min[0] + self.x_offset - self.inking_accel_margin
        self.y_starting_point = self.print_mesh_min[1] + self.y_offset
        self.pass_size = self.print_mesh_max[0] - self.print_mesh_min[0]
        return

# --------------------------------------------------------- #

    def connect_board(self, gcmd):
        gcmd.respond_info("Connecting to printhead...")
        if self.conn:
            gcmd.respond_info("Already connected")
            return
        try:
            self.conn = serial.Serial(self.serial, 115200)
            gcmd.respond_info("Connected!")
        except SerialException:
            gcmd.respond_info("Unable to connect")
        return

    def disconnect_board(self, gcmd):
        gcmd.respond_info("Disconnecting printhead...")
        if self.conn:
            self.conn.close()
            self.conn = None
            gcmd.respond_info("Disconnected")
        return

    def configure_board(self, gcmd):
        gcmd.respond_info("Configuring printhead board...")
        try:
            self.conn.write(f"SET_ACCEL_DELAY {self.inking_accel_delay}\n".encode())
            time.sleep(0.001)
            self.conn.write(f"SET_ROW_DELAY {self.inking_columns_delay}\n".encode())
            time.sleep(0.001)
            gcmd.respond_info("Board configured")
        except SerialException:
            gcmd.respond_info("Communication error")
        return

# --------------------------------------------------------- #

    def get_bounds(self, gcmd):
        gcmd.respond_info("Getting bounds...")
        self.get_2d_print_bounds()
        self.x_starting_point = self.print_mesh_min[0] + self.x_offset - self.inking_accel_margin
        self.y_starting_point = self.print_mesh_min[1] + self.y_offset
        self.pass_size = self.print_mesh_max[0] - self.print_mesh_min[0]
        return

    def get_2d_print_bounds(self):
        exclude_objects = self.printer.lookup_object("exclude_object", None)
        if not exclude_objects:
            return
        objects = exclude_objects.get_status().get("objects", [])
        if not objects:
            return

        xs, ys = [], []
        for obj in objects:
            for point in obj["polygon"]:
                xs.append(point[0])
                ys.append(point[1])

        self.print_mesh_min = [min(xs), min(ys)]
        self.print_mesh_max = [max(xs), max(ys)]
        return
# --------------------------------------------------------- #

    def set_directory(self, gcmd):
        directory = gcmd.get("DIRECTORY")
        if directory:
            self.pictures_directory = self.pictures_base_path + directory
            gcmd.respond_info(f"Directory set to {self.pictures_directory}")
        return

    def set_layer(self, gcmd):
        self.current_layer = gcmd.get_int("LAYER")
        self.load_layer_data(gcmd)
        gcmd.respond_info(f"Layer set to {self.current_layer}")
        return

    def load_layer_data(self, gcmd):
        try:
            with open(f"{self.pictures_directory}/{self.current_layer}.json") as f:
                data = json.load(f)
            self.current_layer_data = data["slices"]
            self.current_layer_passes = len(self.current_layer_data)
            gcmd.respond_info(f"Loaded {self.current_layer_passes} passes")
        except Exception as e:
            gcmd.respond_info(f"Failed to load layer data: {e}")
        return

# --------------------------------------------------------- #

    def write_layer_pass_to_board(self, gcmd, pass_data):
        self.conn.write("CLEAR_PIC\n".encode())
        time.sleep(0.002)
        for i in pass_data:
            #gcmd.respond_info(f"{i}\n")
            self.conn.write(f"SET_PIC {i}\n".encode())
            time.sleep(0.001)
        gcmd.respond_info(f"loaded!")
        self.is_pass_loaded = True
        return

    def trigger_print(self):
        self.gcode.run_script_from_command('SET_PIN PIN=fire_pin VALUE=1')
        self.gcode.run_script_from_command('SET_PIN PIN=fire_pin VALUE=0')
        return

# --------------------------------------------------------- #

    def move_printhead_to_starting_point(self, gcmd):
        pos = self.toolhead.get_position()
        gcmd.respond_info(f"going to z:{self.z_lift}")
        self.toolhead.manual_move([None, None, pos[2] + self.z_lift], self.z_lift_speed)
        self.toolhead.wait_moves()
        gcmd.respond_info(f"going to x:{self.x_starting_point} y:{self.y_starting_point}")
        self.toolhead.manual_move([self.x_starting_point, self.y_starting_point, None], self.positioning_speed)
        return

    def print_layer(self, gcmd):
        self.position_before_2d_printing = self.toolhead.get_position()
        gcmd.respond_info(f"moving to starting point")
        self.move_printhead_to_starting_point(gcmd)
        self.toolhead.wait_moves()

        for ix, pass_data in enumerate(self.current_layer_data):
            self.current_pass = ix
            gcmd.respond_info(f"Loading pass {ix+1}/{self.current_layer_passes}")
            #self.write_layer_pass_to_board(gcmd, pass_data)
            gcmd.respond_info(f"purge")
            #self.move_to_purge(gcmd)
            gcmd.respond_info(f"Printing pass {ix+1}/{self.current_layer_passes}")
            self.print_pass(gcmd)
            if (ix < self.current_layer_passes):
                self.prepare_next_pass(gcmd)

        gcmd.respond_info(f"returning to position_before_2d_printing")
        self.toolhead.manual_move([self.position_before_2d_printing[0], self.position_before_2d_printing[1], None], self.positioning_speed)
        self.toolhead.wait_moves()
        self.toolhead.manual_move([None, None, self.position_before_2d_printing[2]], self.z_lift_speed)
        self.toolhead.wait_moves()
        return

    def print_pass(self, gcmd):
        gcmd.respond_info(f"going to x:{self.x_starting_point + self.pass_size}")
        self.trigger_print()
        self.toolhead.manual_move([self.x_starting_point + self.pass_size, None, None], self.inking_speed)
        self.toolhead.wait_moves()
        return

    def prepare_next_pass(self, gcmd):
        gcmd.respond_info(f"going to x:{self.x_starting_point} y:{self.y_starting_point + (self.inking_pass_height*self.current_pass + 1)}")
        self.toolhead.manual_move([self.x_starting_point, self.y_starting_point + (self.inking_pass_height*(self.current_pass + 1)), None], self.positioning_speed)
        self.toolhead.wait_moves()
        return
    
    def printhead_purge(self, gcmd):
        #self.conn.write(f"PURGE {self.purge_cycles}\n".encode())
        time.sleep(self.purge_delay)
        return
    
    def move_to_purge(self, gcmd):
        pos = self.toolhead.get_position()
        self.toolhead.manual_move([None, None, pos[2] + self.purge_z_lift], self.z_lift_speed)
        self.toolhead.wait_moves()
        self.toolhead.manual_move([self.purge_x_position, self.purge_y_position, None], self.positioning_speed)
        self.toolhead.wait_moves()
        self.printhead_purge(gcmd)
        self.toolhead.manual_move([pos[0], pos[1], pos[2]], self.positioning_speed)
        self.toolhead.wait_moves()
        self.toolhead.manual_move([None, None, pos[2]], self.z_lift_speed)
        self.toolhead.wait_moves()
        return


# --------------------------------------------------------- #

def load_config(config):
    return Inkjet(config)