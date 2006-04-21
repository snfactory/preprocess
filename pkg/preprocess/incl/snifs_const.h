/* === Doxygen Comment ======================================= */
/*! 
 * \file          snifs_const.h
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Fri Apr 21 17:52:47 2006
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \version       0.0
 * \brief         
 *                
 * $Id$
 **/
/* =========================================================== */

#ifdef __cplusplus
extern "C" {
#endif

/* ===== GAINS ================================================== */ 
static const double kGainBlue[2]={0.773,0.744};
//static double kGainRed[2]={0.736,0.792}; // from GL measurments
static const double kGainRed[2]={0.757,0.770};     // YC: correction of gain factor from conts
static const double kGainPhot[4]={1.618,1.576,1.51,1.52};


/* ===== BAD COLUMNS ================================================== */

static const int kBadSectionsBlue=0;
static const int kBadSectionsRedNoVar=6;
static const int kBadSectionsRed=6;
//coordinates are on the final assembled image
static const int kBadSectionsDataRed[kBadSectionsRedNoVar*4]={
  /*  1570,1576,2844,2997,
  1576,1576,1,2843,
  1577,1581,2927,2954,
  1581,1581,1,2926, */
  1017,1018,1,4079,// format is [x1,x2,y1,y2]
  1111,1111,1,3888,
  1029,1029,1,2032,
  1032,1032,1,1000,
  1062,1062,1,620,
  1577,1580,2939,2944,
};
  

  // line 1111 if bad from 1 to 2092 - TBConfirmed
/* very hot
+ 10 pixels to be sure
1570 2915 - 10
1576 2959 
1576 1 2914
WARM
+ 5 pixels here
1577 2937 - 5
1581 2946
1581 1 
Bad column : 
1017 1
1018 4079
*/

static const int kBadSectionsPhot=4;
static const int kBadSectionsDataPhot[kBadSectionsPhot*4]={
  638,638,57,57,
  637,639,58,59,
  640,641,59,59,
  638,640,60,4096};// format is [x1,x2,y1,y2]

/*dead column and spurious stuff at beginning
638 57
637 58 - 639 59
640 59 641 59
638 60 640 4096
*/

/* ===== HOT COLUMNS BLEEDING ============================================ */

static const int kSpecialRed[2]={4,12};
static const double kSpecialCorr[2][12]={{60.7,10.3,8.4,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{61.3,12.1,10.3,7.1,5.1,3.8,3.4,2.6,2.1,1.8,1.2,0.5}};
static const double kSpecialVar[2][12]= {{15.1, 3.1,3.1,1.6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{14.6, 5.4, 5.4,5.4,1.7,0.8,0.8,0.8,0.6,0.5,0.5,0.5}};
static const double kSpecialConservative=2;
static const double kSpecialErrorFast=0.18;
                                         
/* ===== HFFF LINES ================================================== */  

static const int kRedNHfff[2]={11,9};
// starting at 1
  static const int kRedHfffLine[2][11]={
    {  509,1020,1021,2041,2042,2043,2044,2045,2046,3068,3069},
    { 2041,2042,2043,2044,2045,2556,3068,3069,3580,-1,-1}};

static const double kRedHfffVal[2][11]={
  {0.9808,0.9751,1.0138,0.9827,0.9281,0.8741,0.9717,0.9678,0.9925,1.0368,
   0.9918},
  {  0.9852,0.9292,0.8764,0.9695,0.9681,1.0117,1.0267,0.9936,1.0541,1,1}};
static const double kRedHfffSigma=0.004;
  

#ifdef __cplusplus
}
#endif
