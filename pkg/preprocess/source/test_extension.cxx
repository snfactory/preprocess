/* This program to test the functionnalities we need for the multi-extension */

#include "image.hxx"

void test_write(ImageSimple *image, char* keyname, float val){
  char key[lg_name+1];
  image->SetTo(val);
  strcpy(key,keyname);
  image->WrDesc("TESTDESC",CHAR,lg_name+1,key);
  printf("%s *** Written %s, %f\n",image->Name(),keyname,val);
}

void test_read(ImageSimple *image, char* keyname, float val){
  char key[lg_name+1];
  float rval = image->RdFrame(0,0);
  image->RdDesc("TESTDESC",CHAR,lg_name+1,key);
  printf("%s    Expected %s, %f\n",image->Name(),keyname,val);
  printf("%s        read %s, %f\n",image->Name(),key,rval);
 }

void test_create_write_close(ImageSimple *image, char* name, char* keyname, float val){
  image->CreateFrame(name,100,100);
  test_write(image,keyname,val);
  test_read(image,keyname,val);
  image->CloseFrame();
}

void test_create_ext_write_close(ImageSimple *image, char* name, char* keyname, float val,char* extname){
  image->CreateFrame(name,100,100);
  test_write(image,keyname,val);
  test_read(image,keyname,val);
  image->CloseFrame();
}

void test_open_read_close(ImageSimple *image, char* name, char* keyname, float val){
  image->OpenFrame(name);
  test_read(image,keyname,val);
  image->CloseFrame();
}


void test_delete(ImageSimple *image,char * name){
  image->OpenFrame(name,"IO");
  image->DeleteFrame(); 
}


int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -ext1 none -ext2 none");

  // test_extension -in test -ext1 Unnn -ext2 deux
  // test_extension -in test.fits -ext1 Unnn -ext2 deux
  // et avec -noask

  init_session(argv,argc,&arglabel,&argval);

  char baseName[lg_name+1],altBase[lg_name+1],extName1[lg_name+1],extName2[lg_name+1],*name;
  const int nbpattern = 9;
  char pattern[nbpattern][lg_name+1]={"base   ","ReBase  ","Ext1    ",
      "ReExt1  ","ReReBase","ReReExt1","Ext2    ","Ext1(2) ","Ext2(1) "};
  float val[nbpattern]={0.1,1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1};

  strcpy(baseName,argval[0]);
  sprintf(altBase,"%s",argval[0]);
  sprintf(extName1,"%s[%s]",argval[0],argval[1]);
  sprintf(extName2,"%s[%s]",argval[0],argval[2]);
  /* raw constructor */
  ImageSimple *image = new ImageSimple();
  ImageSimple image2;

  /* test the no-extension version */
  /* first create */
  int i=0;
  name=baseName;
  test_create_write_close(image,name,pattern[i],val[i]);
  test_open_read_close(image,altBase,pattern[i],val[i]);
  /* recreate and delete */
  i=1;
  test_create_write_close(image,name,pattern[i],val[i]);
  test_open_read_close(image,altBase,pattern[i],val[i]);
  test_delete(image,name);

  /* create with 1 extension */
  i=2;
  name=extName1;
  test_create_ext_write_close(image,name,pattern[i],val[i],argval[1]);
  test_open_read_close(image,name,pattern[i],val[i]);
  /* recreate */
  i=3;
  test_create_ext_write_close(image,name,pattern[i],val[i],argval[1]);
  test_open_read_close(image,name,pattern[i],val[i]);

  /* recreate base , but with an existing and different file */
  i=4;
  name = baseName;
  test_create_write_close(image,name,pattern[i],val[i]);
  test_open_read_close(image,altBase,pattern[i],val[i]);
  /* recreate ext1 , but with an existing and different file */
  i=5;
  name = extName1;
  test_create_ext_write_close(image,name,pattern[i],val[i],argval[1]);
  image->OpenFrame(name,"IO"); // not possible to open IO and I the same file.
  test_read(image,pattern[i],val[i]);

  /* create a second extension (shall be in 2nd slot)*/
  i=6;
  int k=5;
  name = extName2;
  test_create_ext_write_close(&image2,name,pattern[i],val[i],argval[2]);
  test_read(image,pattern[k],val[k]);
  image->CloseFrame();

  image2.OpenFrame(name,"IO"); // not possible to open IO and I the same file.
  test_read(&image2,pattern[i],val[i]);
  test_open_read_close(image,extName1,pattern[k],val[k]);
  

  /* re-create the first extension while the second extension is opened 
   but the first one already exists */
  i=7;
  k=6;
  name = extName1;
  image->CreateFrame(name,100,100);
  test_write(image,pattern[i],val[i]);
  test_read(&image2,pattern[k],val[k]);
  image2.CloseFrame();
  test_read(image,pattern[i],val[i]);
  
  /* re-create the second extension while the first one is open */
  i=8;
  k=7;
  name = extName2;
  image2.CreateFrame(name,100,100);
  test_write(&image2,pattern[i],val[i]);
  test_read(image,pattern[k],val[k]);
  image->CloseFrame();
  test_read(&image2,pattern[i],val[i]);
  image2.CloseFrame();
  
  /* Delete the first frame only --- and reads the 2nd */
  test_delete(image,extName1);
  test_open_read_close(&image2,extName2,pattern[i],val[i]);
  
  delete image;

  /* shall be enough for the moment */
  exit_session(0);
  
}
