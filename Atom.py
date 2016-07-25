import numpy as np

class Atom():
    def __init__(self,linestring,protein_atom=True):

        ##########legacy
        self.type = None



        ##############









        self.protein_atom = protein_atom
        self.original_linestring = linestring
        self.record_type = None
        self.serialnumber = None
        self.atom_name = None
        self.alt_loc = None
        self.residue_name = None
        self.chain_id = None
        self.residue_number = None
        self.code_res_insertion = None
        self.x_coord = None
        self.y_coord = None
        self.z_coord = None
        self.occupancy = None
        self.tempfac = None
        self.segment_id = None
        self.element_symbol = None
        self.charge = None
        ################
        self.parsestring()

        self.modelsplit = False # if set to true, will write endmodel after this residue and startmodel with a random string
        self.term = False # if set to true, will write a term after this residue when pdb writes out

    def get_coord(self):
        x = float(self.x_coord)
        y = float(self.y_coord)
        z = float(self.z_coord)
        return np.array([x,y,z])




    def set_x_coord(self,value):
        self.x_coord = round(float(value),3)
    def set_y_coord(self,value):
        self.y_coord = round(float(value),3)
    def set_z_coord(self,value):
        self.z_coord = round(float(value),3)


    def parsestring(self): # a datatype key, the linestring, and the columns to look for each type
        linestring = self.original_linestring
        if self.protein_atom == True:
            self.record_type = self.parse(1,linestring,(1,4)) #4
        else:
            self.record_type = self.parse(1,linestring,(1,6)) #6
        self.serialnumber =self.parse(2,linestring,(7,11)) # 5
        self.atom_name =self.parse(3,linestring,(12,16))# 6
        self.alt_loc = self.parse(4,linestring,(17,17))#1
        self.residue_name = self.parse(5,linestring,(18,20))# 3
        self.chain_id = self.parse(6,linestring,(22,22))# 1
        self.res_seq_num = self.parse(7,linestring,(23,26)) # 4
        self.code_res_insertion = self.parse(8,linestring,(27,27)) # 1
        self.x_coord = self.parse(9,linestring,(31,38)) # 8
        self.y_coord = self.parse(10,linestring,(39,46)) # 8
        self.z_coord = self.parse(11,linestring,(47,54)) # 8
        self.occupancy = self.parse(12,linestring,(55,60)) # 6
        self.tempfac = float(self.parse(13,linestring,(61,66))) # 6
        self.segment_id = self.parse(14,linestring,(73,76)) # 4
        self.element_symbol = self.parse(15,linestring,(77,78)) # 2
        self.charge = self.parse(16,linestring,(79,80)) # 2
        ################
        #end standard pdb fields





    def parse(self,datatype,linestring,region):
        value = linestring[region[0]-1:region[1]]
        expected_length = (region[1]+1)-region[0]
        if value.isspace():
            return None
        else:
            value = value.strip()

            padded_value = self.padded(datatype,value,expected_length)
            #change types
            if datatype in [9,10,11,12,13,16]:
                return float(padded_value)
            else:
                return padded_value





    def padleft(self,string,expected_length):
        length = len(string)
        padding = expected_length-length
        for space in range(0,padding):
            string = " "+string
        return string

    def padright(self,string,expected_length):
        length = len(string)
        padding = expected_length-length
        for space in range(0,padding):
            string = string+" "
        return string

    def periodprocess(self,string,expected_length_right,expected_length_left): #"left or right"\

        string = self.checkperiod(string)

        split = string.split(".")

        leftstring = split[0].strip()
        left = self.padleft(leftstring,expected_length_left)

        rightstring = split[1].strip()
        right = self.padright(rightstring,expected_length_right)


        return left+"."+right


    def checkperiod(self,string): #add periods
        if "." in string:
            return string
        else:
            return string+".0"



    def padded(self,datatype,value,expected_length):
        if value == None:
            value = " "
        if value != str:
            value = str(value)
        if datatype == 1:
            if value == "ATOM":
                return "ATOM  "
            elif value == "HETATM" or value == "HETA":
                return "HETATM"
            else:
                print("Choked on datatype "+str(value))
        elif datatype == 2:

            if len(str(value)) > expected_length:
                oldvalue = value
                value = str(hex(int(str(value)))[2:])
                #print("Changing "+str(oldvalue)+" to "+value)

            return self.padleft(value,expected_length)
        elif datatype == 3:
            return self.padleft(value,expected_length)
        elif datatype == 4:
            return self.padleft(value,expected_length)
        elif datatype == 5:
            return self.padleft(value,expected_length)
        elif datatype == 6:
            return self.padleft(value,expected_length)
        elif datatype == 7:
            if len(str(value)) > expected_length:
                oldvalue = value
                value = str(hex(int(str(value))*2)[2:expected_length+3]).upper()
                #print("Changing "+str(oldvalue)+" to "+value)
            return self.padleft(value,expected_length)
        elif datatype == 8:
            return self.padleft(value,expected_length)
        elif datatype == 9:
            return self.padleft(self.periodprocess(value,3,3),expected_length)
        elif datatype == 10:
            return self.padleft(self.periodprocess(value,3,3),expected_length)
        elif datatype == 11:
            return self.padleft(self.periodprocess(value,3,3),expected_length)
        elif datatype == 12:
            return self.padleft(self.periodprocess(value,1,2),expected_length)
        elif datatype == 13:
            return self.padleft(self.periodprocess(value,3,2),expected_length+1) #?? not sure why necessary
        elif datatype == 14:
            return self.padright(value,expected_length)
        elif datatype == 15:
            return self.padleft(value,expected_length)
        elif datatype == 16:
            value = value.strip("-") #charges are not supported
            value = value.strip("+")
            return self.padright(value,expected_length)

    def generate_term(self):
        if self.protein_atom == True:
            return self.padded(1,self.record_type,4)+"  "+"     "+"  "+"    "+self.padded(5,self.residue_name,3)+" "+self.padded(6,self.chain_id,1)+self.padded(7,self.res_seq_num,4)+"\n"



    def generate_entry(self):
        if self.protein_atom == True:
            return self.padded(1,self.record_type,6)+""+self.padded(2,self.serialnumber,5)+""+self.padded(3,self.atom_name,4)+" "+self.padded(5,self.residue_name,3)+" "+self.padded(6,self.chain_id,1)+self.padded(7,self.res_seq_num,4)+self.padded(8,self.code_res_insertion,1)+"   "+self.padded(9,self.x_coord,8)+self.padded(10,self.y_coord,8)+self.padded(11,self.z_coord,8)+self.padded(12,self.occupancy,6)+self.padded(13,self.tempfac,6)+"     "+self.padded(14,self.segment_id,4)+self.padded(15,self.element_symbol,2)+self.padded(16,self.charge,2)+"\n"
        else:
            return self.padded(1,self.record_type,6)+""+self.padded(2,self.serialnumber,5)+""+self.padded(3,self.atom_name,)+" "+self.padded(5,self.residue_name,3)+" "+self.padded(6,self.chain_id,1)+self.padded(7,self.res_seq_num,4)+self.padded(8,self.code_res_insertion,1)+"   "+self.padded(9,self.x_coord,8)+self.padded(10,self.y_coord,8)+self.padded(11,self.z_coord,8)+self.padded(12,self.occupancy,6)+self.padded(13,self.tempfac,6)+"     "+self.padded(14,self.segment_id,4)+self.padded(15,self.element_symbol,2)+self.padded(16,self.charge,2)+"\n"

        # skipped alt_loc entry to accomodate hydrogen double digits

    def print(self,verbose=False):
        if verbose == False:
            print(self.generate_entry().strip("\n"))
        else:
            print("Record Type: "+str(self.record_type))
            print("Atom Number: "+str(self.serialnumber))
            print("Atom Name: "+str(self.atom_name))
            print("Alt Loc: "+str(self.alt_loc))
            print("Residue Name: "+str(self.residue_name))
            print("Chain ID: "+str(self.chain_id))
            print("Residue Seq Number: "+str(self.res_seq_num))
            print("Code Res Insert: "+str(self.code_res_insertion))
            print("X Coordinate: "+str(self.x_coord))
            print("Y Coordinate: "+str(self.y_coord))
            print("Z Coordinate: "+str(self.z_coord))
            print("Occupancy: "+str(self.occupancy))
            print("TempFactor: "+str(self.tempfac))
            print("Segment ID: "+str(self.segment_id))
            print("Element Symbol: "+str(self.element_symbol))
            print("Charge: "+str(self.charge))