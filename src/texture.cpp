#include "texture.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

// uint8_to_float() and float_to_uint8() moved to texture.h so that they can be referenced from software_renderer.cpp

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE(sky): 
  // The starter code allocates the mip levels and generates a level 
  // map simply fills each level with a color that differs from its
  // neighbours'. The reference solution uses trilinear filtering
  // and it will only work when you have mipmaps.

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  for(size_t i = 1; i < tex.mipmap.size(); ++i) {
    MipLevel& currMipLevel = tex.mipmap[i];
    MipLevel& prevMipLevel = tex.mipmap[i-1];
    
    for(size_t x=0; x<currMipLevel.width; ++x)
      for(size_t y=0; y<currMipLevel.height; ++y)
      {
        Color c[2][2], color=Color(0,0,0,0);
        for(size_t dx=0; dx<2; ++dx)
          for(size_t dy=0; dy<2; ++dy)
          {
            size_t px=2*x+dx , py=2*y+dy;
            uint8_to_float( &c[dx][dy].r , &prevMipLevel.texels[4*(py*prevMipLevel.width + px)] );
            color += 0.25*c[dx][dy];
          }
        float_to_uint8( &currMipLevel.texels[4*(y*currMipLevel.width + x)] , &color.r );
      }
  }
}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task ?: Implement nearest neighbour interpolation

  if(level>=(int)tex.mipmap.size())
    return Color(1,0,1,1);
  else
  {
    MipLevel mipLevel=tex.mipmap[level];

    int bx = (int) floor(u * (float) mipLevel.width);
    int by = (int) floor(v * (float) mipLevel.height);

    if(0<=bx and bx<mipLevel.width and 0<=by and by<mipLevel.height)
    {
      Color color;
      uint8_to_float(&color.r,&mipLevel.texels[4 * (by * mipLevel.width + bx)]);
      return color;
    }
    else
      return Color(1,0,1,1);
  }

}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task ?: Implement bilinear filtering

  if(level>=(int)tex.mipmap.size())
    return Color(1,0,1,1);
  else
  {
    MipLevel &mipLevel=tex.mipmap[level];
    if(u<0 or u>mipLevel.width or v<0 or v>mipLevel.height)
      return Color(1,0,1,1);
    
    u = max( u , (float) (1.0/(2.0*mipLevel.width) + 1.0e-15) ); u = min( u , (float) (1.0-1.0/(2.0*mipLevel.width) - 1.0e-15) );
    v = max( v , (float) (1.0/(2.0*mipLevel.height) + 1.0e-15) ); v = min( v , (float) (1.0-1.0/(2.0*mipLevel.height) - 1.0e-15) );

    float x = u * (float) mipLevel.width;
    float y = v * (float) mipLevel.height;
    int fx = (int) floor(x);
    int fy = (int) floor(y);

    int x0,x1;
    if( x-(float)fx > 0.5) {
      x0=fx, x1=fx+1;
    }
    else {
      x0=fx-1; x1=fx;
    }
    int y0,y1;
    if( y-(float)fy > 0.5) {
      y0=fy; y1=fy+1;
    }
    else {
      y0=fy-1; y1=fy;
    }

    Color c[2][2];
    uint8_to_float(&c[0][0].r,&mipLevel.texels[4 * (y0 * mipLevel.width + x0)]);
    uint8_to_float(&c[0][1].r,&mipLevel.texels[4 * (y1 * mipLevel.width + x0)]);
    uint8_to_float(&c[1][0].r,&mipLevel.texels[4 * (y0 * mipLevel.width + x1)]);
    uint8_to_float(&c[1][1].r,&mipLevel.texels[4 * (y1 * mipLevel.width + x1)]);

    Color color = c[0][0] * ( (float)x1 + 0.5 - x ) * ( (float)y1 + 0.5 - y )
                + c[1][0] * ( x - (float)x0 - 0.5 ) * ( (float)y1 + 0.5 - y )
                + c[0][1] * ( (float)x1 + 0.5 - x ) * ( y - (float)y0 - 0.5 )
                + c[1][1] * ( x - (float)x0 - 0.5 ) * ( y - (float)y0 - 0.5 );

    return color;
  }
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 8: Implement trilinear filtering

  float d = log2f(max(u_scale*tex.width,v_scale*tex.height));
  d=max(d,1.0e-15f); d=min(d,(float)(tex.mipmap.size()-1)-1.0e-15f);

  int fd = (int)floor(d);
  return (1.0f+(float)fd-d)*sample_bilinear(tex,u,v,d) + (d-(float)fd)*sample_bilinear(tex,u,v,d+1);
}

} // namespace CMU462
