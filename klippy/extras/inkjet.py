import serial
import json

def load_config(config):
    return Inkjet(config)

class Inkjet:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.reactor = self.printer.get_reactor()
        self.gcode = self.printer.lookup_object('gcode')
        self.toolhead = self.printer.lookup_object('toolhead')

        # --- CONFIGS --- #

        self.serial_port = config.get('serial_port')
        self.base_directory_path = config.get('base_directory_path')
        self.x_offset = config.getfloat('x_offset')
        self.y_offset = config.getfloat("y_offset")
        self.z_offset = config.getfloat("z_offset")
        self.z_lift = config.getfloat("z_lift")
        self.z_lift_speed = config.getfloat("z_lift_speed")
        self.move_speed = config.getfloat("move_speed")
        self.move_accel = config.getfloat("move_accel")
        self.print_speed = config.getfloat("print_speed")
        self.print_accel = config.getfloat("print_accel")
        self.overscan = config.getfloat("overscan")
        self.fire_lead_us = config.getfloat("fire_lead_us")
        self.column_fire_interval_us = config.getfloat("column_fire_interval_us")
        #self.swath_axis = config.get("swath_axis")

        # --- STATES --- #

        self.ser = None
        self.payload_search_directory = None
        self.payload_filename = None
        self.payload_absolute_path = None
        self.payload_fp = None
        self.payload = None
        self.stored_position = None
        self.print_origin = [0.0, 0.0]
        self.pass_stride = 0
        self.swath_length = 0
        self.total_swaths = 0
        self.mesh_min = (0.0, 0.0)
        self.mesh_max = (0.0, 0.0)
        self.current_swath_id = 0
        self.is_printing = False

        self.register_commands()
        self.printer.register_event_handler("klippy:shutdown", self.on_shutdown)

    # --- COMMANDS --- #

    def register_commands(self):
        reg = self.gcode.register_command
        reg("INKJET_CONNECT", self.cmd_connect,desc="")
        reg("INKJET_STATUS", self.cmd_status,desc="")
        reg("INKJET_DISCONNECT", self.cmd_disconnect,desc="")
        reg("INKJET_SET_DIRECTORY", self.cmd_set_directory, desc="")
        reg("INKJET_LOAD", self.cmd_load,desc="")
        reg("INKJET_INFO", self.cmd_info,desc="")
        reg("INKJET_VALIDATE", self.cmd_validate,desc="")
        reg("INKJET_CLOSE", self.cmd_close,desc="")
        reg("INKJET_CLEAR", self.cmd_clear,desc="")
        reg("INKJET_GET_BOUNDS", self.cmd_get_bounds,desc="")
        reg("INKJET_PRINT", self.cmd_print,desc=(""))

    # --- LIFECICLE --- #

    def on_shutdown(self):
        try:
            self.close_serial()
        except Exception:
            pass
        try:
            self.clear()
        except Exception:
            pass
        return

    # --- INKJET_CONNECT --- #

    def open_serial(self):
        try:
            self.ser = serial.Serial(self.serial_port, 115200, timeout=1)
        except Exception as e:
            raise self.gcode.error(f"ERROR")
        return

    def cmd_connect(self, gcmd):
        self.open_serial()
        return

    # --- INKJET_STATUS --- #

    def cmd_status(self, gcmd):
        return

    # --- INKJET_DISCONNECT --- #

    def close_serial(self):
        if self.ser is not None:
            try:
                self.ser.close()
            except Exception as e:
                raise self.gcode.error(f"ERROR")
        return

    def cmd_disconnect(self, gcmd):
        if self.ser is None:
            return
        self.close_serial()
        return

    # --- INKJET_SET_DIRECTORY --- #

    def set_directory(self, directory):
        self.payload_search_directory = directory
        return

    def cmd_set_directory(self, gcmd):
        directory = gcmd.get("DIRECTORY")
        self.set_directory(directory)
        self.gcode.respond_info(f"OK", log=True)
        return

    # --- INKJET_LOAD --- #

    def open_file(self, filename):
        if self.payload_fp is not None:
            self.close_file()
        try:
            self.payload_filename = filename
            self.payload_absolute_path = self.base_directory_path + self.payload_search_directory + self.payload_filename
            self.payload_fp = open(self.payload_absolute_path, "r", encoding="utf-8")
            self.payload = json.load(self.payload_fp)
            self.pass_stride = self.payload["stats"]["pass_stride_mm"]
            self.swath_length = self.payload["stats"]["swath_length_mm"]
            self.total_swaths = self.payload["stats"]["total_swaths"]
        except Exception as e:
            self.payload_fp = None
            raise self.gcode.error(f"ERROR")
        
    def cmd_load(self, gcmd):
        filename = gcmd.get("FILENAME")
        if not filename:
            raise gcmd.error("ERROR")
        self.open_file(filename)
        self.gcode.respond_info(f"OK", log=True)
        return

    # --- INKJET_INFO --- #

    def cmd_info(self, gcmd):
        return

    # --- INKJET_VALIDATE --- #

    def cmd_validate(self, gcmd):
        return

    # --- INKJET_CLOSE --- #

    def close_file(self):
        if self.payload_fp is not None:
            try:
                self.payload_fp.close()
            finally:
                self.payload_fp = None
        self.payload = None

    def cmd_close(self, gcmd):
        self.close_file()
        self.gcode.respond_info(f"OK", log=True)
        return

    # --- INKJET_CLEAR --- #

    def cmd_clear(self, gcmd):
        self.close_file()
        self.payload_filename = None
        self.payload_absolute_path = None
        self.payload = None
        self.stored_position = None
        self.print_origin = [0.0, 0.0]
        self.pass_stride = 0
        self.swath_length = 0
        self.total_swaths = 0
        self.current_swath_id = 0
        self.mesh_min = (0.0, 0.0)
        self.mesh_max = (0.0, 0.0)
        self.is_printing = False
        self.gcode.respond_info(f"OK", log=True)

    # --- INKJET_GET_BOUNDS --- #

    def get_mesh_bounds(self):
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

        self.mesh_min = [min(xs), min(ys)]
        self.mesh_max = [max(xs), max(ys)]
        return

    def get_bounds(self):
        self.get_mesh_bounds()
        x_origin = self.mesh_min[0] + self.x_offset - self.overscan
        y_origin = self.mesh_min[1] + self.y_offset
        self.print_origin = [x_origin, y_origin]

    def cmd_get_bounds(self, gcmd):
        self.get_bounds()
        return

    # --- INKJET_PRINT --- #

    def move_to_swath_start(self,swath_index):
        x_position = self.print_origin[0]
        y_position = self.print_origin[1] + (swath_index * self.pass_stride)
        self.toolhead.manual_move([x_position, y_position, None], self.move_speed)
        self.toolhead.wait_moves()
        return
    
    def move_to_print_origin(self):
        pos = self.toolhead.get_position()
        self.toolhead.manual_move([None, None, pos[2] + self.z_lift], self.z_lift_speed)
        self.toolhead.wait_moves()
        self.move_to_swath_start(0)
        return
    
    def store_position(self):
        self.stored_position = self.toolhead.get_position()
    
    def return_to_stored_position(self):
        self.toolhead.manual_move([self.stored_position[0], self.stored_position[1], None], self.move_speed)
        self.toolhead.wait_moves()
        self.toolhead.manual_move([None, None, self.stored_position[2]], self.z_lift_speed)
        self.toolhead.wait_moves()
        return

    def execute_swath_pass(self):
        self.trigger_print()
        self.toolhead.manual_move([self.print_origin[0] + self.swath_length, None, None], self.print_speed)
        self.toolhead.wait_moves()
        return

    def prepare_next_swath(self):
        self.current_swath_id += 1
        self.load_swath_data()
        self.move_to_swath_start(self.current_swath_id)
        return
    
    def execute_print_job(self):
        self.current_swath_id = 0
        self.store_position()
        self.move_to_print_origin()
        self.load_swath_data()

        for swath_id in range(self.total_swaths):

            self.execute_swath_pass()
            if (swath_id < self.total_swaths - 1):
                self.prepare_next_swath()

        self.return_to_stored_position()
        return
    
    def cmd_print(self, gcmd):
        self.execute_print_job()
        return
    
    # --- CONTROLLER BOARD RELATED --- #

    def trigger_print(self):
        return
    
    def load_swath_data(self):
        return
    
    def set_swath_data(self):
        return
    
# TODO
# Override no print_origin em caso de enviar posicionamentos no INKJET_GET_BOUNDS
# Implementar envio de dados
# Utilizar mesh z (manual_move -> G1)
