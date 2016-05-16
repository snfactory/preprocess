/* === Doxygen Comment ======================================= */
/*! 
 * \file          snifs_defs.h
 * \copyright     (c) 2004 SNFactory collaboration
 * \author        Emmanuel Gangler <e.gangler@ipnl.in2p3.fr>
 * \brief         
 *    The  snifs_defs.h contains some specific definitions and constants 
 *    for SNIFS
 * $Id: snifs_defs.h,v 1.2 2009/11/30 15:04:02 gangler Exp $
 **/
/* =========================================================== */

#ifndef SNIFS_DEFS_H
#define SNIFS_DEFS_H

#include "IFU_io.h"

const int kNChannel = 4;
//const char kChannelName[kNChannel][lg_name+1]={"Blue channel","Red channel","Photometry","Guiding"} ;
const char kChannelName[kNChannel][lg_name+1]={"B","R","P","G"} ;

#endif
