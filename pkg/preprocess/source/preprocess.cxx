
#include "bichip.hxx"
#include "image.hxx"
#include "catorfile.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -bias null -dark null -flat null");
  init_session(argv,argc,&arglabel,&argval);

  char tmp_name[lg_name+1],inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  
  BiChipSnifs * bias=0;
  ImageSnifs *dark=0, *flat=0;
  
  if (is_set(argval[2]))
    bias = new BiChipSnifs(argval[2]);
  if (is_set(argval[3]))
    dark = new ImageSnifs(argval[3]);
  if (is_set(argval[4]))
    flat = new ImageSnifs(argval[4]);

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
  
    BiChipSnifs * in=new BiChipSnifs(inName);
    sprintf(tmp_name,"mem://%s",inName);

    BiChipSnifs * tmp_out= new BiChipSnifs(*in,tmp_name,FLOAT,1);
    ImageSnifs *out = tmp_out->Preprocess(outName,bias,dark,flat);

    delete in;
    delete tmp_out;
    delete out;
  }

  if (flat) delete flat;
  if (dark) delete dark;
  if (bias) delete bias;

  exit_session(0);
  
}
