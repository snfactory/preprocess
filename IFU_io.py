
#global variables
TK = 0
ASK = 0


import items
#/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#!
#!.blk                                   Handling ERRORS
#!
#!.func                     set_control_level()
#!
#!.purp              To define user's control level
#!
#!.desc
#! void set_control_level(level)
#! int level;
#!             FATAL (default value) causes program to exit on error.
#!             WARNING display warning message and return status ( <> 0)
#!             NONE return status ( <> 0)
#!.ed
#-------------------------------------------------------------------- */



Error_Control_Level = [] #   /* error control level */
Error_Current_Level = -1
Erase_File_Sav = 0;     #   /* save of file erase flag */

def set_control_level(level):
    global Error_Control_Level
    global Error_Current_Level
    
    Error_Current_Level += 1
    
    Error_Control_Level.append(level)

#/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#!
#!.func                       Handle_Error()
#!
#!.purp    defines what to do according to set level control.
#!
#!.desc
#! void Handle_Error(routine,status)
#! char *routine;        name of the routine in which error occured
#! int status;           status code
#!.ed
#-------------------------------------------------------------------- */

def Handle_Error(routine,  status):

   
   

    if (Error_Control_Level is not None):
        level = Error_Control_Level[Error_Current_Level];
    
    else:
        level = 'FATAL'

    if (level == 'FATAL'):
                                    # /* What to do when an error occurs ? */
                #/* causes program to exit */

        errtext = f"FATAL error from routine {routine}"
        print(errtext)
        errtext = "" #print_err
        errtext = get_tiger_errmsg(status,errtext)
        print(errtext) #print_err
        # exit_session(status)
        exit(1)

    elif level== 'WARNING' :        #  /* display warning */

        errtext = "WARNING from routine %s" % routine
        print("WARNING: ", errtext) #print warning
        errtext = get_tiger_errmsg(status,errtext)
        print("WARNING: ", errtext)
        return

    elif level== 'NONE' :  #           /* return status */
        return

    else:       #        /* it's up to the user to handle the error */
        return


def disable_user_warnings():
    set_control_level('NONE')




# /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# !
# !.func                         RD_desc()
# !
# !.purp   reads values in descriptor, returns number of values read
# !.desc
# ! nbval = RD_desc(anyfile,descr,descr_type,nb_elt,val)
# !
# ! (type) *anyfile;      image, spectrum or table
# ! char *descr;          name of descriptor (or list of synonym)
# ! short type;           type of data
# ! int nb_elt;           number of elements to read
# ! (type) *val;       array to store values
# !.ed
# ----------------------------------------------------------------------*/
import numpy as np
def RD_desc(anyfile, descr, type,  nb_elt, val) :
  val = np.array(val)

  ltype = type
  nb = nb_elt
  status = 0
  pt_file = items.Anyfile
  item = "" 
  list_of_desc = ""
  errtext = ""

#   list_of_desc = (char *)malloc((strlen(descr)+1)*sizeof(char));
  list_of_desc = descr

  disable_user_warnings()

#   /* Read synomyms up to 1st one present */
  items = list_of_desc.split("|")
  status = Read_one_desc(anyfile, item, ltype, nb, val)
  for item in items:
    if (status < 0):
          break
    if (item is not None):
      status = Read_one_desc(anyfile, item, ltype, nb, val)   
  

  restore_user_warnings()

  if (status < 0):
    pt_file = anyfile #cast to an anyfile?
    print(errtext, "RD_desc: file= %s desc=%s" % (pt_file.name,descr))
    Handle_Error(errtext, 'ERR_BAD_DESC')
  return(status)

def restore_user_warnings():
	Error_Current_Level-=1


def Read_one_desc(anyfile, descr,  type,  nb_elt, val):
  class _MyBreak(Exception): pass
  stat= 0
  nb_val=nb_elt
  nbread = 0
  i = 0
  l = 0
  nulls = 0
  pt_file = None #anyfile
  pt_float = 0.0
  pt_double = 0.00
  errtext = ""
  dsc_items = [] #descr items 
  unit = 0
#ifdef FITS
  fptr = None #fits file
  buffer = ""
  pt_buf = ""
  int_datatype = 0
  n = 0
#endif

  pt_file = anyfile
  try:
        match (pt_file.data_format):
        # TODO TODO TODO find where the functions scdrdc scdrdi .... etc are defined and convert to py
                case 'MIDAS_FORMAT' :
                        match (type):
                                case 'CHAR' :
                                        val = val.astype(str) 
                                        stat = SCDRDC(pt_file.imno,descr,1,1,nb_val,nbread,val,unit, nulls)
                                        if (not stat): 
                                                val = str(val)[0:nbread]
                                case 'SHORT' | 'INT' | 'LONG' :
                                        val = val.astype(int)
                                        stat = SCDRDI(pt_file.imno,descr,1,nb_val, nbread,val,unit, nulls)
                                case 'FLOAT' | 'DOUBLE':
                                        val = val.astype(float)
                                        stat = SCDRDR(pt_file.imno,descr,1,nb_val, nbread,val,unit, nulls)
                                
                case 'IRAF_FORMAT' | 'STSDAS_FORMAT':
                        len_descr = len(descr)
                        # /* check for table descriptors */
                        if (pt_file.file_type != 'T_TABLE'):
                                match (type) :
                                        case 'CHAR' :
                                                uhdgvt(pt_file.imno,descr, one, one, nbread,val,stat,len_descr,nb_val)
                                                
                                        case 'SHORT' :
                                                uhdgvs((pt_file.imno),descr, one, nb_val,nbread,val,stat,len_descr)
                                                
                                        case 'INT' :
                                                uhdgvi( (pt_file.imno),descr, one, nb_val, nbread,val,
                                                stat,len_descr)
                                                
                                        case 'LONG' :
                                                uhdgvl( (pt_file.imno),descr, one, nb_val, nbread,val,
                                                stat,len_descr)
                                                
                                        case 'FLOAT' :
                                                uhdgvr( (pt_file.imno),descr, one, nb_val, nbread,val,
                                                stat,len_descr)
                                                
                                        case 'DOUBLE' :
                                                uhdgvd( (pt_file.imno),descr, one, nb_val, nbread,val,
                                                stat,len_descr)
                                                
                        
                        else:
                                match (type):#		/* table descriptor */
                                        case 'CHAR' :
                                                uthgtt( (pt_file.imno),descr,val, stat,len_descr,nb_val)
                                                
                                        case 'SHORT' | 'LONG' | 'INT' :
                                                uthgti( (pt_file.imno),descr,val, stat,len_descr)
                                                
                                        case 'FLOAT' :
                                                uthgtr( (pt_file.imno),descr,val, stat,len_descr)
                                                
                                        case 'DOUBLE' :
                                                uthgtd( (pt_file.imno),descr,val, stat,len_descr)
                                                
                                

                case 'FITS_A_FORMAT' | 'FITS_B_FORMAT' | 'EURO3D_FORMAT':

                        fptr = pt_file.external_info
                        stat = 0
                        nbread = nb_val
                        if (nbread == 0):
                                raise _MyBreak

                        if (type == 'CHAR'):
                                val = ''
                                int_datatype = get_datatype_code(pt_file.data_format,type)
                                if  (descr == "COMMENT" or  descr == "HISTORY"):
                                        fits_read_key_longstr(fptr, descr,  pt_buf, None,  stat)
                                        if (stat): raise _MyBreak
                                        n = min(len(pt_buf),nb_val)
                                        val = pt_buf[:n]
                                else:
                                        fits_read_key(fptr, int_datatype, descr, buffer, None,  stat)
                                        if (stat): raise _MyBreak
                                        if (strlen(buffer) < nb_val):
                                                val = buffer
                                        else:
                                                val=buffer[:nb_val]
                                raise _MyBreak
                        if (nb_val == 1):
                                int_datatype = get_datatype_code(pt_file.data_format,type)
                                fits_read_key(fptr, int_datatype, descr, val, None,  stat)
                                raise _MyBreak
                
                        match (type):
                                case 'CHAR' :
                                        fits_read_keys_str(fptr,descr,1,nb_val,val, nbread, stat)
                                case 'SHORT' | 'LONG' | 'INT' :
                                        fits_read_keys_lng(fptr,descr,1,nb_val,val, nbread, stat)
                                case 'FLOAT' :
                                        fits_read_keys_flt(fptr,descr,1,nb_val,val, nbread, stat)
                                case 'DOUBLE' :
                                        fits_read_keys_dbl(fptr,descr,1,nb_val,val, nbread, stat)
                case 'TIGER_FORMAT' :
                        stat = 0
                        if (pt_file.external_info == None):
                                stat = -1
                                raise _MyBreak
                        dsc_items = pt_file.external_info
                        i=0
                        while (i<dsc_items.nb_descr and dsc_items.descr_list[i].descr_name!=descr):
                                i += 1
                        if (i == dsc_items.nb_descr):
                                stat = 'ERR_NODESC'
                                raise _MyBreak
                
                        if (dsc_items.descr_list[i].data_type == 'LONG' and (type == 'INT')):
                                dsc_items.descr_list[i].data_type = type
                        if ((dsc_items.descr_list[i].data_type == 'INT') and (type == 'LONG')):
                                dsc_items.descr_list[i].data_type = type

                        nb_val = min(nb_val,dsc_items.descr_list[i].nb_values)
                        nbread = nb_val

                        if (dsc_items.descr_list[i].data_type != type):
                                if ((dsc_items.descr_list[i].data_type == 'DOUBLE') and (type == 'FLOAT')):
                                        pt_float = val.astype(float)
                                        for l in range(nb_val):
                                                pt_float[l] = float(dsc_items.descr_list[i].descr_value.d_data[l])
                        
                                elif ((dsc_items.descr_list[i].data_type == 'FLOAT') and (type == 'DOUBLE')):
                                        pt_double = val.astype(float)
                                        for l in range(nb_val):
                                                pt_double[l] = float(dsc_items.descr_list[i].descr_value.f_data[l])
                                else:
                                        stat = 'ERR_BAD_DESC'
        
         
                        else:
                                val = dsc_items.descr_list[i].descr_value.c_data[:nb_val]
                        
  except _MyBreak():
        pass
  if (stat):
        print(errtext, "RD_desc: file= %s desc=%s"% (pt_file.name,descr))
        stat = get_tiger_errcode(pt_file.data_format,stat)
        Handle_Error(errtext, stat)
  if (stat < 0):
        return(stat)
  else:
        return(nbread)


#/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#!
#!.blk                  Routines for CATALOG read/write
#!
#!.func                         RD_catalog()
#!
#!.purp             successive reads of catalog members
#!.desc
#! RD_catalog(catalog,filename)
#!
#! char *catalog;        name of catalog
#! char *filename;       name of catalog member
#!.ed
#----------------------------------------------------------------------*/

def RD_catalog(catalog, filename):

        # char text[132];
        text = ""
        fd_catal = [] #list of file ptr?  File**
        nameof_catal = [] #list of string
        nb_catal = 0 #static int

        i = 0
        status = 0
        filename = '' #no accessible null terminator in py

        if (nb_catal != 0):
            while (i<nb_catal and  (catalog != nameof_catal[i])):
                    i += 1
        # else:
        #     fd_catal = (FILE **)malloc(sizeof(FILE *));
        #     nameof_catal = (char **)malloc(sizeof(char *));
            
        
        if (i >= nb_catal):
            #fd_catal= (FILE **)realloc((char *)fd_catal,(nb_catal+1)*sizeof(FILE *));
            fd_catal = [None] * (nb_catal + 1)
            nameof_catal = [""] * (nb_catal + 1) 
            #nameof_catal= (char **)realloc((char *)nameof_catal,(nb_catal+1)*sizeof(char *));
            #dont need malloc
            #nameof_catal[nb_catal] = (char *)malloc(lg_name*sizeof(char));
            
            nameof_catal[nb_catal] = catalog
            
            fd_catal[nb_catal] = open(catalog,"r") #TODO handle IO
            
            if (fd_catal[nb_catal] is None): #ie == NULL
                text=f"RD_catalog {catalog}"
                Handle_Error(text,'ERR_BAD_CAT')
                return False
            
            else:
                text = fd_catal[i].readLine() #/* get catalog description line */
                nb_catal+=1
        
        if (text == ''):
                fd_catal[i].close()
        else:
                text = fd_catal[i].readLine()
                if (text == ''):
                    fd_catal[i].close()
                else:
                    # text[first_blk(text)] = '\0'
                    filename = text.strip().split()[0]
                
        
        if (filename == ''):
            nb_catal = 0
            i =0
            return "" #falsy
        
        return filename #truthy



#/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#!
#!.func                   get_tiger_errcode()
#!
#!.purp        returns the error code according to Tiger conventions
#!
#!.desc
#! int get_tiger_errcode(data_format,stat)
#! short data_format;    data format
#! int stat;             error status code
#!
#!.ed
#-------------------------------------------------------------------- */
def get_tiger_errcode(data_format, stat):

    match data_format:
        case 'FITS_A_FORMAT' | 'FITS_B_FORMAT' | 'EURO3D_FORMAT':

    #case FITS_A_FORMAT :
    #case FITS_B_FORMAT :
    #case EURO3D_FORMAT :
            match stat:
                case  101 | 103 | 104:
                    return('ERR_OPEN')
                case  105:
                        return('ERR_CREAT')
                case  106:
                        return('ERR_WRIT')
                case  107 | 108:
                        return('ERR_READ')
                case  110:
                        return('ERR_CLOSE')
                case  111:
                        return('ERR_ALLOC')
                case  112:
                        return('ERR_ACCESS')
                case  201:
                        return('ERR_BAD_HEAD')
                case  202:
                        return('ERR_NODESC')
                case  203:
                        return('ERR_HEADER_SIZE')
                case  204 | 205| 206 | 207 | 208|209|210|211|212|213|214|215|216|217|218 :
            
                        return('ERR_BAD_DESC')
                case  219:
                        return('ERR_NOCOL')
                case  220:
                        return('ERR_BAD_DESC')
                case  221|222|223|224|225|226|227|228|229|230|231|232:
                        return('ERR_BAD_HEAD')
                case  233:
                        return('ERR_IMA_HEAD')
                case  234:
                        return('ERR_BAD_DESC')
                case  235:
                        return('ERR_TBL_HEAD')
                case  236:
                        return('ERR_HEADER_SIZE')
                case  237:
                        return('ERR_NOCOL')
                case  241:
                        return('ERR_BAD_DESC')
                case  251|252|253|261|262|263|301:
                        return('ERR_BAD_HEAD')
                case  302:
                        return('ERR_COL_NUM')
                case  304 | 306:
                        return('ERR_NODATA')
                case  307|308:
                        return('ERR_BAD_PARAM')
                case  309|310|311|312:
                        return('ERR_BAD_TYPE')
                case  314:
                        return('ERR_NODATA')
                case  317|320|321:
                        return('ERR_BAD_PARAM')
                case  322|323:
                        return('ERR_BAD_DESC')
                case  401|402|403|404|405|406|407|408|409|410:
                        return('ERR_BAD_TYPE')
                case  411:
                        return('ERR_BAD_PARAM')
                        return('ERR_OVERFLOW')
                case  501|502|503:
                        return('ERR_BAD_PARAM')
                case  505:
                        return('ERR_BAD_DESC')
                case _: return(stat)
        

        case 'MIDAS_FORMAT':
                match (stat):

                        case -4 : return('ERR_NODATA');
                        case -6 : return('ERR_NODATA')
                        case  0 : return('OK')
                        case  1 : return('ERR_NODESC')
                        case  2 : return('ERR_OVERFLOW')
                        case  6 : return('ERR_ACCESS')
                        case  7 : return('ERR_BAD_PARAM')
                        case  8 : return('ERR_OVERFLOW')
                        case  9 : return('ERR_BAD_DESC')
                        case 15 : return('ERR_BAD_CAT')
                        case 20 : return('ERR_OVERFLOW')
                        case 21 : return('ERR_ALLOC')
                        case 22 : return('ERR_ALLOC')
                        case 23 : return('ERR_OVERFLOW')
                        case 24 : return('ERR_NOTBL')
                        case 25 : return('ERR_COL_NUM')
                        case 26 : return('ERR_ROW_NUM') 
                        case 27 : return('ERR_NOIDENT')
                        case 28 : return('ERR_BAD_COL')
                        case 29 : return('ERR_NOIMPL')
                        case 31 : return('ERR_REN_TBL')
                        case 32 : return('ERR_NOCOL')
                        case _: return(stat)
                
        case 'STSDAS_FORMAT'|'IRAF_FORMAT' :
            match (stat) :
                case -2 : return('ERR_EOF')
                case  0 : return('OK')
                case  6 : return('ERR_BAD_PARAM')
                case  7 : return('ERR_NB_PARAM')
                case 10 : return('ERR_OPEN')
                case 11 : return('ERR_CREAT')
                case 12 : return('ERR_CREAT')
                case 13 : return('ERR_NAXIS')
                case 14 : return('ERR_NAXIS')
                case 15 : return('ERR_BAD_TYPE')
                case 16 : return('ERR_OFFSET')
                case 17 : return('ERR_ACCESS')
                case 18 : return('ERR_CLOSE')
                case 19 : return('ERR_BAD_IMA')
                case 20 : return('ERR_IMA_BOUND')
                case 21 : return('ERR_IMA_EXT')
                case 22 : return('ERR_IMA_BOUND')
                case 23 : return('ERR_READ')
                case 24 : return('ERR_WRIT')
                case 25 : return('ERR_NODESC')
                case 26 : return('ERR_NODESC')
                case 40 : return('ERR_NODESC')
                case 41 : return('ERR_BAD_HEAD')
                case 42 : return('ERR_HEADER_SIZE')
                case 43 : return('ERR_BAD_PARAM')
                case 46 : return('ERR_BAD_PARAM')
                case 47 : return('ERR_NB_PARAM')
                case 48 : return('ERR_DEL_DESC')
                case 49 : return('ERR_DEL_DESC')
                case 50 : return('ERR_NOIMA')
                case 51 : return('ERR_DEL_IMA')
                case 52 : return('ERR_REN_IMA')
                case 53 : return('ERR_IMA_HEAD')
                case 54 : return('ERR_BAD_DESC')
                case 70 : return('ERR_GRAPH_DEV')
                case 101 : return('ERR_ALLOC')
                case 102 : return('ERR_FREE')
                case 112 : return('ERR_BAD_SIZE')
                case _: return(stat)
            
        case 'TIGER_FORMAT' :
            return(stat)
                
     
    return('ERR_FORMAT')

# /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# !
# !.func                   get_tiger_errmsg()
# !
# !.purp     returns the error message corresponding to the given code
# !
# !.desc
# ! void get_tiger_errmsg(stat,msg)
# ! int stat;             error status code
# ! char *msg;            corresponding error message
# !
# !.ed
# -------------------------------------------------------------------- */
def get_tiger_errmsg(stat, msg):

    if(not TK):
        match(stat):
            case 'OK' : msg="Successful return\n"
            case 'UNKNOWN' : msg=''
            case 'ERR_EOF' : msg = "EOF reached\n"
            case 'ERR_BAD_PARAM' : msg = "Input element is bad\n"
            case 'ERR_NB_PARAM' : msg = "Bad input number of elements\n"
            case 'ERR_BAD_TYPE' : msg = "Invalid data type\n"
            case 'ERR_BAD_SIZE' : msg = "Field width not wide enough\n"
            case 'ERR_OPEN' : msg = "Error opening file\n"
            case 'ERR_CREAT' : msg = "Error opening new file\n"
            case 'ERR_READ' : msg = "Error reading file\n"
            case 'ERR_WRIT' : msg = "Error writing to file\n"
            case 'ERR_CLOSE' : msg = "Error closing file\n"
            case 'ERR_ACCESS' : msg = "Invalid image access mode\n"
            case 'ERR_NAXIS' : msg = "Invalid NAXIS parameter\n"
            case 'ERR_NOIDENT' : msg = "Identifier not found\n"
            case 'ERR_OFFSET' : msg = "Error returning offset to data\n"
            case 'ERR_NODATA' : msg = "No data available\n"
            case 'ERR_BAD_IMA' : msg = "Not an image !\n"
            case 'ERR_IMA_BOUND' : msg = "Bad section specification for image\n"
            case 'ERR_IMA_EXT' : msg = "Bad extension for image\n"
            case 'ERR_NOIMA' : msg = "Image does not exist\n"
            case 'ERR_DEL_IMA' : msg = "Error deleting image\n"
            case 'ERR_REN_IMA' : msg = "Error renaming image\n"
            case 'ERR_IMA_HEAD' : msg = "Illegal image header\n"
            case 'ERR_BAD_TBL' : msg = "Not a table !\n"
            case 'ERR_TBL_EXT' : msg = "Bad extension for table\n"
            case 'ERR_NOTBL' : msg = "Table does not exist\n"
            case 'ERR_DEL_TBL' : msg = "Error deleting table\n"
            case 'ERR_REN_TBL' : msg = "Error renaming table\n"
            case 'ERR_TBL_HEAD' : msg = "Illegal table header\n"
            case 'ERR_BAD_COL' : msg = "Error in column format\n"
            case 'ERR_NOCOL' : msg = "Column does not exist\n"
            case 'ERR_COL_NUM' : msg = "Wrong column number\n"
            case 'ERR_ROW_NUM' : msg = "Wrong row number\n"
            case 'ERR_NODESC' : msg = "Header parameter not found\n"
            case 'ERR_BAD_HEAD' : msg = "Illegal data type for header parameter\n"
            case 'ERR_DEL_DESC' : msg = "Cannot delete descriptor\n"
            case 'ERR_BAD_DESC' : msg = "Descriptor bad\n"
            case 'ERR_HEADER_SIZE' : msg = "Out of space in header\n"
            case 'ERR_BAD_CAT' : msg = "Not a catalog !\n"
            case 'ERR_GRAPH_DEV' : msg = "Bad graphics device\n"
            case 'ERR_ALLOC' : msg = "Error allocating dynamic memory\n"
            case 'ERR_FREE' : msg = "Error freeing dynamic memory\n"
            case 'ERR_NOIMPL' : msg = "Not yet implemented\n"
            case 'ERR_OVERFLOW' : msg = "Overflow (column/frame)\n"
            case 'ERR_FORMAT' : msg = "Unknown data format\n"
            case _ : msg=""
        
    else:
        match(stat):
            case 'OK' : msg = "Successful return"
            case 'UNKNOWN' : msg=''
            case 'ERR_EOF' : msg = "EOF reached"
            case 'ERR_BAD_PARAM' : msg = "Input element is bad"
            case 'ERR_NB_PARAM' : msg = "Bad input number of elements"
            case 'ERR_BAD_TYPE' : msg = "Invalid data type"
            case 'ERR_BAD_SIZE' : msg = "Field width not wide enough"
            case 'ERR_OPEN' : msg = "Error opening file"
            case 'ERR_CREAT' : msg = "Error opening new file"
            case 'ERR_READ' : msg = "Error reading file"
            case 'ERR_WRIT' : msg = "Error writing to file"
            case 'ERR_CLOSE' : msg = "Error closing file"
            case 'ERR_ACCESS' : msg = "Invalid image access mode"
            case 'ERR_NAXIS' : msg = "Invalid NAXIS parameter"
            case 'ERR_NOIDENT' : msg = "Identifier not found"
            case 'ERR_OFFSET' : msg = "Error returning offset to data"
            case 'ERR_NODATA' : msg = "No data available"
            case 'ERR_BAD_IMA' : msg = "Not an image !"
            case 'ERR_IMA_BOUND' :  msg = "Bad section specification for image"
            case 'ERR_IMA_EXT' : msg = "Bad extension for image"
            case 'ERR_NOIMA' : msg = "Image does not exist"
            case 'ERR_DEL_IMA' : msg = "Error deleting image"
            case 'ERR_REN_IMA' : msg = "Error renaming image"
            case 'ERR_IMA_HEAD' : msg = "Illegal image header"
            case 'ERR_BAD_TBL' : msg = "Not a table !"
            case 'ERR_TBL_EXT' : msg = "Bad extension for table"
            case 'ERR_NOTBL' : msg = "Table does not exist"
            case 'ERR_DEL_TBL' : msg = "Error deleting table"
            case 'ERR_REN_TBL' : msg = "Error renaming table"
            case 'ERR_TBL_HEAD' : msg = "Illegal table header"
            case 'ERR_BAD_COL' : msg = "Error in column format"
            case 'ERR_NOCOL' : msg = "Column does not exist"
            case 'ERR_COL_NUM' : msg = "Wrong column number"
            case 'ERR_ROW_NUM' : msg = "Wrong row number"
            case 'ERR_NODESC' : msg = "Header parameter not found"
            case 'ERR_BAD_HEAD' :
                            msg = "Illegal data type for header parameter"
            case 'ERR_DEL_DESC' : msg = "Cannot delete descriptor"
            case 'ERR_BAD_DESC' : msg = "Descriptor bad"
            case 'ERR_HEADER_SIZE' : msg = "Out of space in header"
            case 'ERR_BAD_CAT' : msg = "Not a catalog !"
            case 'ERR_GRAPH_DEV' : msg = "Bad graphics device"
            case 'ERR_ALLOC' : msg = "Error allocating dynamic memory"
            case 'ERR_FREE' : msg = "Error freeing dynamic memory"
            case 'ERR_NOIMPL' : msg = "Not yet implemented"
            case 'ERR_OVERFLOW' : msg = "Overflow (column/frame)"
            case 'ERR_FORMAT' : msg = "Unknown data format\n"
            case _: msg = ''
    return msg
        