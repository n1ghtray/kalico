import serial
import logging

from serial import SerialException

class Ink:
    def __init__(self, config):
        self.printer = config.get_printer()
        self.reactor = self.printer.get_reactor()
        self.gcode = self.printer.lookup_object("gcode")

        #printhead related
        self.x_offset = config.getfloat("x_offset")
        self.y_offset = config.getfloat("y_offset")
        self.lift_z = config.getfloat("lift_z")
        self.pass_size = config.getfloat("pass_size")
        self.accel_distance = config.getfloat("accel_distance")
        self.accel_delay = config.get("accel_delay")
        self.rows_delay = config.get("rows_delay")

        #2d printing related
        self.is_printing = False
        self.printing_axis = config.get("printing_axis")

        #serial connection related
        self.con = None
        self.serial_port = config.get("serial")
        self.trigger_pin = config.get("trigger_pin")

        #2d printing data related
        self.pictures_base_path = config.get("pictures_base_path") #/home/pi/printer_data/gcode_ink/__GCODE_FILENAME__


        self.gcode.register_command("INK_CONNECT", self.Connect)
        self.gcode.register_command("INK_DISCONNECT", self.Disconnect)


    def Connect(self, gcmd):
        if self.serial:
            gcmd.respond_info("Already connected!")
            return
        try:
            self.serial = serial.Serial(self.serial_port, 115200) #timeout=0, write_timeout=0
        except SerialException:
            gcmd.respond_info("Unable to connect!")
            return
        
    def Disconnect(self, gcmd=None):
        self.gcode.respond_info("Disconnecting!")
        if self.serial:
            self.serial.close()
            self.serial = None



    #variáveis que precisamos: um trigger quando está tudo ready
    #um trigger quando um gcode é carregado (pode ser um gcode dentro do arquivo .gcode)
    #um trigger quando há mudança de layer (pode ser um gcode dentro do arquivo .gcode)
    #o resto é de boa, teoricamente
    #métodos para implementar:
    #carregar arquivo json da respectiva camada
    #carregar para a ram do placa os respectivos dados de cada passada
    #posicionar e mover a printhead em cada passada
    #método para quando crasha
    #método para quando starta, enviar os dados de ocnfig para a placa