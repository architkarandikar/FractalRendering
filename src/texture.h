#ifndef CMU462_TEXTURE_H
#define CMU462_TEXTURE_H

#include "CMU462.h"

using namespace std;

namespace CMU462 {

static const int kMaxMipLevels = 14;

typedef enum SampleMethod{
  NEAREST,
  BILINEAR,
  TRILINEAR
} SampleMethod;

struct MipLevel {
  size_t width; 
  size_t height;
  std::vector<unsigned char> texels;
};

struct Texture {
  size_t width;
  size_t height;
  std::vector<MipLevel> mipmap;
};

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

class Sampler2D {
 public:

  Sampler2D( SampleMethod method ) : method ( method ) { }

  ~Sampler2D();

  virtual void generate_mips( Texture& tex, int startLevel ) = 0;

  virtual Color sample_nearest(Texture& tex, 
                               float u, float v, 
                               int level = 0) = 0;

  virtual Color sample_bilinear(Texture& tex, 
                                float u, float v, 
                                int level = 0) = 0;

  virtual Color sample_trilinear(Texture& tex, 
                                 float u, float v, 
                                 float u_scale, float v_scale) = 0;
  
  inline SampleMethod get_sample_method() const {
    return method;
  }
 
 protected:

  SampleMethod method;

}; // class Sampler2D

class Sampler2DImp : public Sampler2D {
 public:

  Sampler2DImp( SampleMethod method = TRILINEAR ) : Sampler2D ( method ) { }
  
  void generate_mips( Texture& tex, int startLevel );

  Color sample_nearest(Texture& tex, 
                       float u, float v, 
                       int level = 0);

  Color sample_bilinear(Texture& tex, 
                        float u, float v, 
                        int level = 0);

  Color sample_trilinear(Texture& tex, 
                         float u, float v, 
                         float u_scale, float v_scale);
  
}; // class sampler2DImp

class Sampler2DRef : public Sampler2D {
 public:

  Sampler2DRef( SampleMethod method = TRILINEAR ) : Sampler2D ( method ) { }
  
  void generate_mips( Texture& tex, int startLevel );

  Color sample_nearest(Texture& tex, 
                       float u, float v, 
                       int level = 0);

  Color sample_bilinear(Texture& tex, 
                        float u, float v, 
                        int level = 0);

  Color sample_trilinear(Texture& tex, 
                         float u, float v, 
                         float u_scale, float v_scale);
  
}; // class sampler2DRef

} // namespace CMU462

#endif // CMU462_TEXTURE_H
