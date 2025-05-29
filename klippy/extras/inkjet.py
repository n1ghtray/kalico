# This file may be distributed under the terms of the GNU GPLv3 license.

import serial
import logging
import json
import time

from serial import SerialException

'''
INKJET_CONNECT
INKJET_ENABLE
INKJET_SET_DIRECTORY=dir
INKJET_SET_LAYER=0
INKJET_PRINT
INKJET_PURGE
INKJET_DISABLE
INKJET_DISCONNECT
'''


'''
How does this works?
A 


'''




class Inkjet:
    def __init__(self, config):
                



        self.printer = config.get_printer()
        self.reactor = self.printer.get_reactor()

        self.printer.register_event_handler("klippy:ready", self._handle_connect)

        self.gcode = self.printer.lookup_object("gcode")

        self.toolhead = None
        self.enabled = False

        #3d printing related
        self.current_layer = None


        #printhead related
        self.x_offset = config.getfloat("x_offset", 10) #x offset between printhead and toolhead
        self.y_offset = config.getfloat("y_offset", 10) #y offset between printhead and toolhead
        self.z_offset = config.getfloat("z_offset", 10) #z offset between printhead and 3d print

        self.z_lift = config.getfloat("z_lift", 1) #z lift between non-inking movements
        self.z_lift_speed = config.getfloat("z_lift_speed", 20) #z lift between non-inking movements


        self.inking_speed = config.getfloat("inking_speed", 50) #inking speed
        self.inking_accel = config.getfloat("inking_accel", 50) #inking acceleration
        self.inking_accel_margin = config.getfloat("inking_accel_margin", 10) #margin to printhead accelerate
        self.inking_axis = config.get("inking_axis", "x") #axis for printhead passes
        self.inking_pass_height = config.getfloat("inking_pass_height", 5.4) #pass height

        #printhead controller board
        self.inking_accel_delay = config.getfloat("inking_accel_delay", 100) #time that it takes to achive constant speed
        self.inking_columns_delay = config.getfloat("inking_columns_delay", 100) #delay between each columns of each pass

        self.inking_is_printing = False
        self.inking_start_layer = config.getint("inking_start_layer", 2) #the first layer that enables inking

        self.positioning_speed = config.getfloat("positioning_speed", 100) #speed for positioning movements
        self.positioning_accel = config.getfloat("positioning_accel", 1000) #acceleration for positioning movements
        self.position_before_2d_printing = None

        self.conn = None
        self.serial = config.get("serial", "test") #serial port path
        self.trigger_pin = config.get("trigger_pin", "test")

        self.pictures_base_path = config.get("pictures_base_path", "home/pi/printer_data/gcode_ink/")
        self.pictures_directory = None
        self.current_layer_data = None
        self.current_layer_passes = None
        self.current_pass = None
        self.is_pass_loaded = False

        self.print_mesh_min = None
        self.print_mesh_max = None

        self.gcode.register_command("INKJET_CONNECT", self.connect)
        self.gcode.register_command("INKJET_ENABLE", self.enable)
        self.gcode.register_command("INKJET_CONFIGURE", self.configure)
        self.gcode.register_command("INKJET_DISABLE", self.disable)
        self.gcode.register_command("INKJET_DISCONNECT", self.disconnect)
        self.gcode.register_command("INKJET_SET_DIRECTORY", self.set_directory)
        self.gcode.register_command("INKJET_SET_LAYER", self.set_layer)
        self.gcode.register_command("INKJET_PRINT", self.print_layer)

    def _handle_connect(self):
        self.toolhead = self.printer.lookup_object("toolhead")

    def connect(self, gcmd):
        gcmd.respond_info("connect")
        if self.conn:
            gcmd.respond_info("Already connected")
            return
        try:
            self.conn = serial.Serial(self.serial, 115200) #timeout=0, write_timeout=0
            gcmd.respond_info("Connected!")
        except SerialException:
            gcmd.respond_info("Unable to connect")
            return
        
    def enable(self, gcmd):
        gcmd.respond_info("enable")
        self.enabled = True
        self.get_2d_print_bounds(gcmd)
        return
        
    def configure(self, gcmd):
        gcmd.respond_info("configure")
        #write settings to the printhead controller board
        try:
            self.conn.write("SET_ACCEL_DELAY {}".format(self.accel_delay).encode())
            time.sleep(0.001)
            self.conn.write("SET_ROW_DELAY {}".format(self.rows_delay).encode())
            time.sleep(0.001)
            gcmd.respond_info("Board configurated")
        except SerialException:
            gcmd.respond_info("Unable to communicate")
        return
    
    def disable(self, gcmd):
        gcmd.respond_info("disable")
        self.enabled = False
        return

    def disconnect(self, gcmd):
        gcmd.respond_info("disconnect")
        if self.conn:
            self.conn.close()
            self.conn = None
            gcmd.respond_info("Disconnected")
        return

    def set_directory(self, gcmd):
        gcmd.respond_info("set_directory")
        directory = gcmd.get("DIRECTORY")
        if (directory):
            self.pictures_directory = self.pictures_base_path + directory
        return

    def set_layer(self, gcmd):
        gcmd.respond_info("set_layer")
        layer = gcmd.get_int("LAYER")
        self.current_layer = layer
        self.load_layer_data()
        return

    def load_layer_data(self, gcmd):
        gcmd.respond_info("load_layer_data")
        #load the layer json data
        with open("{}/{}.json".format(self.pictures_directory, str(self.current_layer))) as f:
            data = json.load(f)
        
        self.current_layer_passes = len(data["slices"])
        self.current_layer_data = data["slices"]
        return

    def write_layer_pass_to_board(self, gcmd, pass_data):
        gcmd.respond_info("write_layer_pass_to_board")
        #send the pass data to the printhead controller board
        self.conn.write("CLEAR_PIC\n".encode())
        time.sleep(0.001)
        for i in pass_data:
            self.conn.write("SET_PIC 0x{} \n".format(i).encode())
            time.sleep(0.001)
        self.is_pass_loaded = True
        return
    
    def move_printhead_to_starting_point(self,gcmd):
        gcmd.respond_info("move_printhead_to_starting_point")
        #position the printhead at the starting position for the first pass
        first_pass_starting_x_coord = self.print_mesh_min[0] + self.x_offset - self.accel_margin
        first_pass_starting_y_coord = self.print_mesh_min[1] + self.y_offset
        self.toolhead.manual_move([None, None, self.lift_z], self.z_lift_speed)
        self.toolhead.manual_move([first_pass_starting_x_coord, first_pass_starting_y_coord, None], self.positioning_speed)
        return
    
    def trigger_print(self,gcmd):
        gcmd.respond_info("trigger_print")
        #set the PIN connected to the printhead controller board at HIGH, and then set it LOW
        return

    def print_layer(self,gcmd):
        gcmd.respond_info("print_layer")
        #print an entire layer
        self.position_before_2d_printing = self.toolhead.get_position()
        self.move_printhead_to_starting_point(gcmd)
        for pass_data,ix in self.current_layer_data:
             gcmd.respond_info("loading pass %s" % (ix))
             self.write_layer_pass_to_board(pass_data)
             self.print_pass()
             if ix != len(self.current_layer_data):
                self.prepare_next_pass()
        self.toolhead.manual_move([self.position_before_2d_printing[0], self.position_before_2d_printing[1], None], self.positioning_speed)
        self.toolhead.manual_move([None, None, self.position_before_2d_printing[2]], self.z_lift_speed)

    def purge(self,gcmd):
        return

    def print_pass(self,gcmd):
        gcmd.respond_info("print_pass")
        self.trigger_print()
        printhead_starting_position = self.toolhead.get_position()
        printhead_pass_end_position = [printhead_starting_position[0] + (self.print_mesh_max[0] - self.print_mesh_min[0]),None,None]
        self.toolhead.manual_move([printhead_pass_end_position], self.inking_speed)
        return
    
    def prepare_next_pass(self,gcmd):
        gcmd.respond_info("prepare_next_pass")
        return

    def get_2d_print_bounds(self, gcmd):
        gcmd.respond_info("get_2d_print_bounds")
        exclude_objects = self.printer.lookup_object("exclude_object", None)
        if exclude_objects is None:
            return False
        objects = exclude_objects.get_status().get("objects", [])
        if not objects:
            return False
        # List all exclude_object points by axis and iterate over
        # all polygon points, and pick the min and max or each axis
        list_of_xs = []
        list_of_ys = []
        gcmd.respond_info("Found %s objects" % (len(objects)))
        for obj in objects:
            for point in obj["polygon"]:
                list_of_xs.append(point[0])
                list_of_ys.append(point[1])
        # Define bounds of adaptive mesh area
        self.print_mesh_min = [min(list_of_xs), min(list_of_ys)]
        self.print_mesh_max = [max(list_of_xs), max(list_of_ys)]
        gcmd.respond_info("mesh_min: %s mesh_max: %s" % (self.print_mesh_min, self.print_mesh_max))
        return

def load_config(config):
    return Inkjet(config)


#para testes, vamos fazer o seguinte: um retângulo impresso, com dimensões  77mm x 55mm x 50mm, começando a imprimir 2d na camada 4, altura de camada 0.1mm.
#teoricamente, a imagem precisa ter, sem considerar as margens devido às 3 cores separadas, 909 pixels de largura (300 DPI) e 650 (300 DPI) pixels de altura.
#teremos 496 imagens, uma para cada camada impressa.


#return_pos = self.toolhead.get_position()
#self.toolhead.manual_move([x, y, None], speed)