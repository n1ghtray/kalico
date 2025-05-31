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
        self.inking_accel_delay = config.getfloat("inking_accel_delay", 100000)
        self.inking_columns_delay = config.getfloat("inking_columns_delay", 64)
        self.inking_start_layer = config.getint("inking_start_layer", 0)
        self.positioning_speed = config.getfloat("positioning_speed", 100)
        self.positioning_accel = config.getfloat("positioning_accel", 1000)

        # State and data
        self.is_inking = False
        self.position_before_inking = None
        self.start_pos = None
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
        self.gcode.register_command("INKJET_PURGE", self.purge)
        self.gcode.register_command("INKJET_SET_DIRECTORY", self.set_directory)
        self.gcode.register_command("INKJET_SET_LAYER", self.set_layer)
        self.gcode.register_command("INKJET_PRINT", self.print_layer)
        self.gcode.register_command("FILL_DEBUG_INFO", self.fill_debug_info)

    def _handle_connect(self):
        self.toolhead = self.printer.lookup_object("toolhead")

    def fill_debug_info(self, gcmd):
        self.connect_board(gcmd)
        self.configure_board(gcmd)
        self.pictures_directory = "/home/pi/printer_data/gcode_ink/TEST"
        self.current_layer = gcmd.get_int("LAYER", 300)
        self.load_layer_data(gcmd)
        xmin = gcmd.get_int("XMIN", 10)
        xmax = gcmd.get_int("XMAX", 10)        
        ymin = gcmd.get_int("YMIN", 87)
        ymax = gcmd.get_int("YMAX", 65)
        self.print_mesh_min = [xmin,ymin]
        self.print_mesh_max = [xmax,ymax]
        return

# --------------------------------------------------------- #

    def connect_board(self, gcmd):
        gcmd.respond_info("Connecting to printhead...")
        if self.conn:
            gcmd.respond_info("Already connected")
            return

    def disconnect_board(self, gcmd):
        gcmd.respond_info("Disconnecting printhead...")
        if self.conn:
            self.conn.close()
            self.conn = None
            gcmd.respond_info("Disconnected")
        try:
            self.conn = serial.Serial(self.serial, 115200)
            gcmd.respond_info("Connected!")
        except SerialException:
            gcmd.respond_info("Unable to connect")

    def configure_board(self, gcmd):
        gcmd.respond_info("Configuring printhead board...")
        try:
            self.conn.write(f"SET_ACCEL_DELAY {self.inking_accel_delay}\n".encode())
            time.sleep(0.001)
            self.conn.write(f"SET_ROW_DELAY {self.inking_columns_delay}\n".encode())
            gcmd.respond_info("Board configured")
        except SerialException:
            gcmd.respond_info("Communication error")

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

# --------------------------------------------------------- #

    def set_directory(self, gcmd):
        directory = gcmd.get("DIRECTORY")
        if directory:
            self.pictures_directory = self.pictures_base_path + directory
            gcmd.respond_info(f"Directory set to {self.pictures_directory}")

    def set_layer(self, gcmd):
        self.current_layer = gcmd.get_int("LAYER")
        self.load_layer_data(gcmd)
        gcmd.respond_info(f"Layer set to {self.current_layer}")

    def load_layer_data(self, gcmd):
        try:
            with open(f"{self.pictures_directory}/{self.current_layer}.json") as f:
                data = json.load(f)
            self.current_layer_data = data["slices"]
            self.current_layer_passes = len(self.current_layer_data)
            gcmd.respond_info(f"Loaded {self.current_layer_passes} passes")
        except Exception as e:
            gcmd.respond_info(f"Failed to load layer data: {e}")

# --------------------------------------------------------- #

    def write_layer_pass_to_board(self, gcmd, pass_data):
        self.conn.write("CLEAR_PIC\n".encode())
        time.sleep(0.001)
        for i in pass_data:
            self.conn.write(f"SET_PIC 0x{i}\n".encode())
            time.sleep(0.001)
        self.is_pass_loaded = True

    def trigger_print(self):
        self.gcode.run_script_from_command('SET_PIN PIN=fire_pin VALUE=1')
        self.gcode.run_script_from_command('SET_PIN PIN=fire_pin VALUE=0')
        pass

# --------------------------------------------------------- #

    def move_printhead_to_starting_point(self):
        self.toolhead.manual_move([None, None, self.z_lift], self.z_lift_speed)
        self.toolhead.move([self.x_starting_point, self.y_starting_point, None], self.positioning_speed)

    def print_layer(self, gcmd):
        self.position_before_2d_printing = self.toolhead.get_position()
        self.move_printhead_to_starting_point()

        for ix, pass_data in enumerate(self.current_layer_data):
            self.current_pass = ix
            gcmd.respond_info(f"Loading pass {ix+1}/{self.current_layer_passes}")
            self.write_layer_pass_to_board(gcmd, pass_data)
            gcmd.respond_info(f"Printing pass {ix+1}/{self.current_layer_passes}")

            self.print_pass()
            if ix + 1 < self.current_layer_passes:
                self.prepare_next_pass()

        self.toolhead.move([self.position_before_2d_printing[0], self.position_before_2d_printing[1], None], self.positioning_speed)
        self.toolhead.move([None, None, self.position_before_2d_printing[2]], self.z_lift_speed)

    def print_pass(self):
        self.trigger_print()
        self.toolhead.manual_move([self.pass_size, None, None], self.inking_speed)

    def prepare_next_pass(self):
        self.toolhead.manual_move([-self.start_pos[0], self.inking_pass_height*self.current_pass, None], self.positioning_speed)
        return

    def purge(self):
        return

# --------------------------------------------------------- #

def load_config(config):
    return Inkjet(config)