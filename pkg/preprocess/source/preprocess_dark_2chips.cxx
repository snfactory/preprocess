
#include "bichip.hxx"
#include "imagesnifs.hxx"
#include "catorfile.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -bias null");
  init_session(argv,argc,&arglabel,&argval);

  char tmp_name[lg_name+1],inName[lg_name+1],outName[lg_name+1];

  CatOrFile inCat(argval[0]);
  CatOrFile outCat(argval[1]);
  
  BiChipSnifs * bias=0;
  
  // Load once auxilliary files
  if (is_set(argval[2]))
    bias = new BiChipSnifs(argval[2]);

  while (inCat.NextFile(inName) && outCat.NextFile(outName)) {
  
    BiChipSnifs * in=new BiChipSnifs(inName);

    BiChipSnifs * tmp_out= new BiChipSnifs(*in,outName,FLOAT,1);
    delete in;
  
    tmp_out->PreprocessBias();
    if (bias)
      tmp_out->SubstractBias();
    tmp_out->AddPoissonNoise()

    delete tmp_out;
  }

  if (bias) delete bias;

  exit_session(0);
  
}
