/* === Doxygen Comment ======================================= */
/*! 
 * \file          image.cxx
 * \copyright     (c) 2003 CRAL-Observatoire de Lyon
 * \date          Wed Aug  6 17:44:35 2003
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#include "ioslice.hxx"
#include "fitsio.h"
#include "longnam.h"
#include "data_io.h"
extern IO_Format OutputIO;
/*
Logic for the alloc :
1) open/create the image
2) allocate memory for the requested MaxLines
   fMaxLines is the unmber of AllocatedLines
3) load lines (triggered by Wr or Rd Frame)
*/

/* performances :
copy_image : 1056*4128 -> idem  : 2.2 CPU s.
copy_image_slice : -nlines 4128 : 2.2 CPU s.
copy_image_slice : -nlines    1 : 2.2 CPU s.
*/

/* ##### IO SLICE ################################################# */

/* ===== constructor/destructor ======================================= */

/* ----- IoSlice -------------------------------------------------- */
IoSlice::IoSlice(int MaxLines) : IoPlain() {
  fAllocLines=0;
  fFirstLine=0;
  fLastLine=0;
  SetMaxLines(MaxLines);
}

/* ----- ~ImageSimple ------------------------------------------------- */
IoSlice::~IoSlice(){
  if (fLoaded)
    CloseFrame();
}

/* ===== Overwritten virtual methods ================================== */

/* ----- WrFrame ------------------------------------------------------ */
void IoSlice::WrFrame(int col, int line, double value){
  //Wrapper to WR_frame
  if (line<FirstLine() || line>=LastLine())
    LoadFromLine(line,0);
  WR_frame(Frame(),col,line-FirstLine(),value);
}

/* ----- RdFrame ------------------------------------------------------ */
double IoSlice::RdFrame(int col, int line) {
  //Wrapper to RD_frame
  if (line<FirstLine() || line>=LastLine())
    LoadFromLine(line,1);
  return (double) RD_frame(Frame(),col,line-FirstLine());
}

/* ----- OpenFrame ---------------------------------------------------- */
int IoSlice::OpenFrame(char *name, char *mode){
  // Wrapper to open_frame
  if (fLoaded)
    print_error("IoSlice::OpenFrame Try to reopen a frame\n");
  fLoaded=true;
  int ret = (!header_frame(Frame(),name,mode) && MemAlloc());
  fExists = new int[Ny()];
  for (int i=0;i<Ny();i++)
    fExists[i]=1;
  return ret;
}

/* ----- CloseFrame --------------------------------------------------- */
int IoSlice::CloseFrame(){
  // Wrapper to close_frame
  WriteCurrent();
  MemFree();
  if (fLoaded) {
    fLoaded=false;
    delete[] fExists;
    return close_frame(Frame());
  } else {
    return 0;
  }
}

/* ----- DeleteFrame --------------------------------------------------- */
int IoSlice::DeleteFrame(){
  // Wrapper to delete_frame
  // Note that the file should be closed first. (at least in iolibs 6.0b++)
  if (fLoaded)
    print_error("IoSlice::DeleteFrame image %s is not closed ...",Name());
  return delete_frame(Frame());
}


/* ----- CreateFrame -------------------------------------------------- */
int IoSlice::CreateFrame(char *Name,int Nx, int Ny, short Type ){
  // this CreateFrame is less powerfull than teh original, but
  // this is because I don't understand all parameters !
  int npix[2];
  double start[2];
  double step[2];

  start[0]=npix[0]=Nx;
  start[1]=npix[1]=Ny;
  step[0]=step[1]=1;
  fLoaded=true;

  int ret=create_frame(Frame(),Name,npix,start,step,Type,"","");
  free_frame_mem(Frame());
  MemAlloc();
  fExists = new int[Ny];
  for (int i=0;i<Ny;i++)
    fExists[i]=0;
  return ret;
}

/* ----- CreateFrameFull -------------------------------------------------- */
int IoSlice::CreateFrameFull(char *name,int *Npix,double* Start, double* Step, short Type,char* Ident, char* U) {

  fLoaded=true;
  int ret = create_frame(Frame(),name,Npix,Start,Step,Type,Ident,U);
  free_frame_mem(Frame());
  MemAlloc();
  fExists = new int[Ny()];
  for (int i=0;i<Ny();i++)
    fExists[i]=0;

  return ret;
  
}


/* ===== Specific methods            ================================== */

/* ----- AllocMem -------------------------------------------------- */
int IoSlice::MemAlloc(){
  if (!fLoaded)
    print_error("IoSlice::MemAlloc image is not yet loaded");
  if (fMaxLines==kIoAll)
    fMaxLines = Ny();
  if (fMaxLines) {
    int ny = Frame()->ny;
    if (Frame()->ny > fMaxLines)
      Frame()->ny = fMaxLines;
    int ret = alloc_frame_mem(Frame(),Frame()->data_type);
    fAllocLines = Frame()->ny;
    Frame()->ny = ny;
    // initialise line reading
    fFirstLine = 0; 
    fLastLine = 0; 
    return ret;
  }
  return OK;
}

/* ----- MemFree -------------------------------------------------- */
int IoSlice::MemFree(){
  if (!fLoaded)
    print_error("IoSlice::MemFree cannot free a non-loaded image");
  if (fAllocLines) {
    fAllocLines=0;
    fFirstLine=0;
    fLastLine=0;
    return free_frame_mem(Frame());
  }
  
  return OK;
}

/* ----- MemIsAlloc -------------------------------------------------- */
int IoSlice::MemIsAlloc(){
  return (fAllocLines >0);
}

/* ----- SetMaxLines ---------------------------------------------- */
void IoSlice::SetMaxLines(int Max){
  // does the alloc if the image is loaded
  fMaxLines = Max;
  if (fLoaded && fAllocLines != (fMaxLines > Ny() ? Ny() :fMaxLines ) ) {
    if (fAllocLines) {
      if (fLastLine>fFirstLine)
        WriteCurrent();
      MemFree();
    }
    
    MemAlloc();
  }
}


/* ----- LoadFromLine ---------------------------------------------- */
void IoSlice::LoadFromLine(int Line,int Read){
  // triggered by the RdFrame and the WrFrame
  if (fLastLine>fFirstLine)
    WriteCurrent();
  if (Line < 0 || Line >=Ny()) {
    print_msg("IoSlice::LoadFromLine requested line %d not available",Line);
    fFirstLine=fLastLine=0;
  }

  // all tests passed - do it.
  fFirstLine = Line;
  if (Ny()>fFirstLine+fAllocLines)
    fLastLine = fFirstLine+fAllocLines ;
  else
    fLastLine = Ny();
  // load only existing lines
  int first=-111; // arbitrary value to placate compiler
  for (int line = fFirstLine; line < fLastLine; line++ ) {
    while (line <fLastLine && !fExists[line])
      line++;
    first = line;
    while (line < fLastLine && fExists[line])
      line++;
    if (first<fLastLine)
      LoadLineRange(first,line);
  }
  if (Read && first != fFirstLine)
    print_warning("IoSlice::LoadFromLine : you may try to read nonexsting lines");
}

/* ----- LoadFromLine ---------------------------------------------- */
void IoSlice::LoadLineRange(int First, int Last){
  // low-level routine to load the lines

  long group=0L;
  long npix = Nx()*(Last-First);
  long firstpix = Nx()*First+1L;
  int nbread,status=0;
  int offset = (First-fFirstLine) * Nx ();
  fitsfile *fptr = (fitsfile*)Frame()->external_info;
  switch (Frame()->data_type) {
    case SHORT :
      if (fits_read_img_sht(fptr,group,firstpix,npix,(short)0,
                            Frame()->data.s_data + offset,&nbread,&status)) {
        status = ERR_READ; break;
      }
      break;
    case LONG :
    case INT :
      if (fits_read_img_lng(fptr,group,firstpix,npix,(int)0,
                          Frame()->data.l_data + offset,&nbread,&status)) {
        status = ERR_READ; break;
      }
      break;
    case FLOAT :
      if (fits_read_img_flt(fptr,group,firstpix,npix,(float)0,
                            Frame()->data.f_data + offset,&nbread,&status)) {
        status = ERR_READ; break;
      }
      break;
    case DOUBLE :
      if (fits_read_img_dbl(fptr,group,firstpix,npix,(double)0,
                            Frame()->data.d_data + offset,&nbread,&status)) {
        status = ERR_READ; break;
      }
      break;
  }
  if (status) {
    char errtext[lg_name+1];
    sprintf(errtext,"IoSlice::LoadLineRange frame %s lines %d-%d",Name(),First,Last);
    status = get_tiger_errcode(Frame()->data_format,status);
    Handle_Error(errtext,status);
  }
}

/* ----- WriteCurrent ---------------------------------------------- */
void IoSlice::WriteCurrent(){

  if (Frame()->iomode == (int)I_MODE)
    return;
  if (fFirstLine>=fLastLine)
    return;

  int datatype = get_datatype_code(OutputIO.basic_io,Frame()->data_type);
  long npix = Nx()*(fLastLine-fFirstLine);
  long firstpix = Nx()*fFirstLine + 1L;
  int status=0;
  fitsfile * fptr = (fitsfile *) Frame()->external_info;
  if (fits_write_img(fptr,datatype,firstpix,npix,
                     Frame()->data.s_data,&status)) {
    status = ERR_WRIT;
  }
  if (status) {
    char errtext[lg_name+1];
    sprintf(errtext,"close_frame: frame %s",Name());
    status = get_tiger_errcode(Frame()->data_format,status);
    Handle_Error(errtext,status);
  }
  for (int i=fFirstLine;i<fLastLine;i++)
    fExists[i]=1;
  
}
