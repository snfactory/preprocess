/* === Doxygen Comment ======================================= */
/*! 
 * \file          snifs_const.h
 * \copyright     (c) 2003 SNIFS-Supernova Factory Experiment
 * \date          Tue Jul  5 14:13:13 2005
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

static const double kGainBlue[2]={0.773,0.744};
//static double kGainRed[2]={0.736,0.792}; // from GL measurments
static const double kGainRed[2]={0.757,0.770};     // YC: correction of gain factor from conts
static const double kGainPhot[4]={1.618,1.576,1.51,1.52};


static const int kBadSectionsBlue=0;
static const int kBadSectionsRed=5;
//coordinates are on the final assembled image
static const int kBadSectionsDataRed[kBadSectionsRed*4]={
  1570,1576,2844,2997,
  1576,1576,1,2843,
  1577,1581,2927,2954,
  1581,1581,1,2926,
  1017,1018,1,4079};// format is [x1,x2,y1,y2]
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

#ifdef __cplusplus
}
#endif
