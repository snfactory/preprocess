
#include "imagesnifs.hxx"
#include "catorfile.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  int xdir=1,ydir=1;
  double zscale=1.0;

  set_arglist("-in none -out none -xflip -yflip -zscale null");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  CatOrFile catOut(argval[1]);
  if (is_true(argval[2]))
    xdir = -1;
  if (is_true(argval[3]))
    ydir = -1;
  if (is_set(argval[4]))
    get_argval(5,"%lf",&zscale);

  char inName[lg_name+1],outName[lg_name+1];
  while (catIn.NextFile(inName) && catOut.NextFile(outName)) {
    print_msg("Opening %s",inName);
    ImageSnifs *in = new ImageSnifs(inName);
    ImageSnifs *out = new ImageSnifs(*in,outName,0,0);
              
    out->ImportFlip(in,xdir,ydir,zscale);
    
    delete in;
    delete out;
  }
  
  exit_session(0);
  
}
