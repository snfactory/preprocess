#include "image.hxx"
#include "catorfile.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -key none -val none -type INT|FLOAT|CHAR");
  init_session(argv,argc,&arglabel,&argval);

  char storage[lg_name+1];
  int type,length=1;
  if (!strcasecmp(argval[3],"int")) {
    get_argval(2,"%d",storage);
    type=INT;
  }
  if (!strcasecmp(argval[3],"float")) {
    get_argval(2,"%f",storage);
    type=FLOAT;
  }
  if (!strcasecmp(argval[3],"char")) {
    strcpy(storage,argval[2]);
    type=CHAR;
    length = lg_name+1;
  }


  char inName[lg_name+1];  
  CatOrFile catIn(argval[0]);
  while (catIn.NextFile(inName)){
    ImageSimple *in = new ImageSimple(argval[0],"IO");
    in->WrDesc(argval[1],type,length,storage);
    delete in;
  }
  
  exit_session(0);
}
