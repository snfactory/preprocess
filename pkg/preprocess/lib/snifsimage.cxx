
#include "snifsimage.h"

SnifsImage SnifsImage::SnifsImage () {
  image = new IMAGE2D;
}

SnifsImage::~SnifsImage() {
  delete image;  
}
