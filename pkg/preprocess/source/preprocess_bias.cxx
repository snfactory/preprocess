
#include "bichip.hxx"
#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -stackout null -sigcut 3.0 -nlines 5000");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  CatOrFile catOut(argval[1]);
  int nlines;
  get_argval(4,"%d",&nlines);

  if (is_set(argval[2])) {
    BiChipStackSnifs * in=new BiChipStackSnifs(&catIn,"I",nlines);// nlines = 1
    // has to be benchmarked ...
    // in that case BiChipStackSnifs::PreprocessBias may have to expand
    // the image with a SetNLines before the copy.
    BiChipStackSnifs * tmp= in->PreprocessBias(&catOut,nlines);
    delete in;

    double sigma;
    get_argval(3,"%lf", &sigma);
    KGauss k(sigma);
    BiChipSnifs* out = tmp->Kombine(argval[2],&k);
    delete tmp;
    delete out;

  } else { // no need to store temporary information
    char inName[lg_name+1],outName[lg_name+1];
    while (catIn.NextFile(inName) && catOut.NextFile(outName)) {
      print_msg("Opening %s",inName);
      BiChipSnifs * in=new BiChipSnifs(inName);
      BiChipSnifs * out= new BiChipSnifs(*in,outName,FLOAT,1);
      delete in;
      out->PreprocessBias();

      delete out;  
    }
  }

  exit_session(0);
}
