#include "image.hxx"
#include "analyser.hxx"
#include "catorfile.hxx"
#include "section.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1];
  double sigcut = 5.0;

  set_arglist("-in none -sec [4:1027,1:4102] -sigcut 5.0");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  Section *sec = new Section(argval[1]);
  sscanf(argval[2],"%f",&sigcut);

  printf("Name Mean Sigma PixOut PixOutMean Max\n");

  while(inCat.NextFile(inName)) {
    ImageSnifs in(inName);
    ImageAnalyser ana(&in,sec);
    double min,max;

    in.MinMax(sec,&min,&max);
    
    printf("%s %f %f %d %f %f \n",in.Name(), 
           ana.MeanLevel(), sqrt(ana.StatsVariance()),ana.NPixOut(sigcut), ana.OutPixMean(sigcut),max);
  }
  delete sec;

  exit_session(0);  
}
