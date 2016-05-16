/* === Doxygen Comment ======================================= */
/*! 
 * \file          fclass_snifs.h
 * \copyright     (c) 2004 SNFactory collaboration
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \brief         
 *    The  fclass_snifs.h contains the specific fclasses for the
 *    SNIFS instrument reduction and analysis software
 *    that are not yet included in the IFU_C_iolibs or any other
 *    distributed 3D LCL library
 * $Id: fclass_snifs.h,v 1.7 2015/03/10 11:13:43 ycopin Exp $
 **/
/* =========================================================== */

/* Renaming of some unclear values in gendef.h */

#define	RAW_DARK_FRAME	DARK_FRAME	/* Raw dark frame */
//#define RAW_DOME_FRAME	RAW_DOM_FRAME	/* Raw Dome Frame */
//#define PRE_DOME_FRAME	RAW_DOM_FRAME	/* Pre Dome Frame */

/* 5xx numbers are for reduction up to the calibrated data cube */

#define PRE_DARK_FRAME	70	/* Preprocessed dark frame */

#define RAW_SCA_FRAME	60 /* Raw SCALA frame (see snifs_keywords.h) */
#define PRE_SCA_FRAME	61 /* Preprocessed SCALA frame (see snifs_keywords.h) */

/* numbers underneath may be obsolete now ... (E.G. 05.27.10)*/
/* 6xx numbers are for datacube extracted informations */

#define FCLASS_TO_BE_IMPLEMENTED 600

/* 9xx are for specific photometric row or processed frame */

#define RAW_PHOTO_FOCUS  901  /* output for a photometric focus run */

