#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 4 (part 2): 
  // Set svg to normalized device coordinate transformation. Your input
  // arguments are defined as SVG canvans coordinates.
  this->x = x;
  this->y = y;
  this->span = span; 

  Matrix3x3 m;
  m(0,0) = 0.5/span; m(0,1) = 0.0; m(0,2) = -0.5*x/span + 0.5;
  m(1,0) = 0.0; m(1,1) = 0.5/span; m(1,2) = -0.5*y/span + 0.5;
  m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0;
  set_canvas_to_norm(m);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CMU462