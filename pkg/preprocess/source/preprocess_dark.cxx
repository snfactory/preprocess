
//#include "bichip.hxx"
#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"
#include "preprocessor.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  int nlines;
  
  set_arglist("-in none -out none -stackout null -bias null -sigma 4.0 -nlines 5000 -all");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  CatOrFile catOut(argval[1]);

  ImageSnifs * bias=0;  
  // Load once auxilliary files
  if (is_set(argval[3]))
    bias = new ImageSnifs(argval[3]);

  get_argval(5,"%d",&nlines);
  
  Preprocessor P;
  if (is_true(argval[6]))
    P.SetAllImage(1);

  char inName[lg_name+1],outName[lg_name+1];

  if (is_set(argval[2])) {
    ImageStackSnifs * tmp = new ImageStackSnifs(nlines);

    P.SetIoMethod(kIoSlice);
    while (catIn.NextFile(inName) && catOut.NextFile(outName)) {
      print_msg("Opening %s",inName);
      ImageSnifs * outFile = P.PreprocessDark(inName,outName,bias);
      tmp->AddImage(outFile);
    }

    double sigma;
    get_argval(4,"%lf", &sigma);
    Kombinator * k = new KGaussPoisson(sigma);
    ImageSnifs *out = tmp->Kombine(argval[2],k);
    delete k;
    delete out;
    delete tmp;

  } else { // no need to store temporary information

    while (catIn.NextFile(inName) && catOut.NextFile(outName)) {
      print_msg("Opening %s",inName);
      ImageSnifs * out = P.PreprocessDark(inName,outName,bias);
      delete out;  
    }
  }
  
  if (bias) 
    delete bias;

  exit_session(0);
}
