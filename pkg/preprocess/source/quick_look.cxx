#include "image.hxx"
#include "analyser.hxx"
#include "catorfile.hxx"
#include "section.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  char inName[lg_name+1];
  double sigcut, limit;

  set_arglist("-in none -sec [4:1027,1:4102] -sigcut 5.0 -limit 10000");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile inCat(argval[0]);
  Section *sec = new Section(argval[1]);
  sscanf(argval[2],"%lf",&sigcut);
  sscanf(argval[3],"%lf",&limit);

  printf("Name Mean Sigma PixOut PixOutMean PixOver Max\n");

  while(inCat.NextFile(inName)) {
    ImageSnifs in(inName);
    char tmp_name[lg_name+1];
    sprintf(tmp_name,"mem://%s",argval[1]);
    ImageSnifs out(inName,tmp_name,FLOAT,1);
    out.SubstractOverscan();
    ImageAnalyser ana(&out,sec);
    double min,max;

    out.MinMax(sec,&min,&max);
    
    printf("%s %f %f %d %f %d %f \n",in.Name(), 
           ana.MeanLevel(), sqrt(ana.StatsVariance()),ana.NPixOut(sigcut), ana.OutPixMean(sigcut),ana.NPixOver(limit), max);
  }
  delete sec;

  exit_session(0);  
}
