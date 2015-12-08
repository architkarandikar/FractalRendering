#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {

inline void alpha_compost( unsigned char* dst, Color newColor ) {
  Color prevColor;
  uint8_to_float( &prevColor.r , dst );

  Color color;
  color.a = 1.0f - (1.0f-newColor.a) * (1.0f-prevColor.a);
  color.r = (1.0f-newColor.a) * prevColor.r + newColor.r ;
  color.g = (1.0f-newColor.a) * prevColor.g + newColor.g ;
  color.b = (1.0f-newColor.a) * prevColor.b + newColor.b ;

  float_to_uint8( dst , &color.r );
}

// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  supersample_target = new unsigned char[4 * sample_rate * sample_rate * target_w * target_h];
  memset(supersample_target, 255, 4 * sample_rate * sample_rate * target_w * target_h);

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    group_transformation = canvas_to_screen;
    draw_element(svg.elements[i]);
  }

  // set top level transformation
  transformation = canvas_to_screen;

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y++;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y++;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y--;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y--;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

  delete[] supersample_target;
}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 3: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;

  // draw_svg() has been modified for supersampling support
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 3: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

  // draw_svg() has been modified for supersampling support
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 4 (part 1):
  // Modify this to implement the transformation stack

  transformation = group_transformation * (element->transform);

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      group_transformation = group_transformation * (element->transform);
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

  delete[] supersample_target;
}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit

  //cerr<<ellipse.center.x<<" "<<ellipse.center.y<<" :: "<<ellipse.radius.x<<" "<<ellipse.radius.y<<"\n";
  //cerr<<ellipse.style.fillColor<<"\n";

  rasterize_ellipse( ellipse.center.x, ellipse.center.y, ellipse.radius.x, ellipse.radius.y, ellipse.style.fillColor );
}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  for( int dy = 0; dy < sample_rate; ++dy )
    for( int dx = 0; dx < sample_rate; ++dx ) {
      alpha_compost( &supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] , color );
    }
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 1: 
  // Implement line rasterization

  float eps = (float) 1.0e-15;

  if( abs(x0-x1) < eps && abs(y0-y1) < eps )
    rasterize_point( x0, y0, color );
  else if( abs(x0-x1) > abs(y0-y1) ) {
    if( x0 > x1 ) {
      swap(x0,x1);
      swap(y0,y1);
    }
    for( int x = (int) floor(x0); x < (int) floor(x1); ++x ) {
      float y = y0 + ( ( (float) x - x0 ) * ( y1 - y0 ) / ( x1 - x0 ) );
      rasterize_point( (float) x, y, color );
    }
  }
  else {
    if( y0 > y1 ) {
      swap(x0,x1);
      swap(y0,y1);
    }
    for( int y = (int) floor(y0); y < (int) floor(y1); ++y ) {
      float x = x0 + ( ( (float) y - y0 ) * ( x1 - x0 ) / ( y1 - y0 ) );
      rasterize_point( x, (float) y, color );
    }
  }
}

float SoftwareRendererImp::line_side_test( float x0, float y0,
                                          float x1, float y1,
                                          float x, float y ) {
  return ( ( x - x0 ) * ( y1 - y0 ) - ( y - y0 ) * ( x1 - x0 ) );
}


void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 2: 
  // Implement triangle rasterization

  float eps = (float) 1.0e-15;

  int bxmin = min( min( (int) floor(x0), (int) floor(x1) ), (int) floor(x2)  );
  int bxmax = max( max( (int) ceil(x0), (int) ceil(x1) ), (int) ceil(x2)  );
  int bymin = min( min( (int) floor(y0), (int) floor(y1) ), (int) floor(y2)  );
  int bymax = max( max( (int) ceil(y0), (int) ceil(y1) ), (int) ceil(y2)  );

  for (int y = max(bymin,0); y <= min(bymax,(int)target_h-1); ++y )
    for( int x = max(bxmin,0); x <= min(bxmax,(int)target_w-1); ++x )
      for( int dy = 0; dy < sample_rate; ++dy )
        for( int dx = 0; dx < sample_rate; ++dx ) {
          float cx = (float) x + 1.0 * (2*dx+1) / (2*sample_rate), cy = (float) y + 1.0 * (2*dy+1) / (2*sample_rate);
          if( line_side_test( x0, y0, x1, y1, cx, cy ) * line_side_test( x0, y0, x1, y1, x2, y2 ) > -eps && 
              line_side_test( x1, y1, x2, y2, cx, cy ) * line_side_test( x1, y1, x2, y2, x0, y0 ) > -eps && 
              line_side_test( x2, y2, x0, y0, cx, cy ) * line_side_test( x2, y2, x0, y0, x1, y1 ) > -eps ) {
            alpha_compost( &supersample_target[4 * ( (x + y * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] , color );
          }
        }
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task ?: 
  // Implement image rasterization

  for(int y = max((int)ceil(y0),0); y<=min((int)floor(y1),(int)target_h-1); ++y)
    for(int x = max((int)ceil(x0),0); x<=min((int)floor(x1),(int)target_w-1); ++x)
      for( int dy = 0; dy < sample_rate; ++dy )
        for( int dx = 0; dx < sample_rate; ++dx )
        {
          float cx = (float) x + 1.0 * (2*dx+1) / (2*sample_rate), cy = (float) y + 1.0 * (2*dy+1) / (2*sample_rate);
          float tx = (cx-x0) / (x1-x0), ty = (cy-y0) / (y1-y0);
          Color color = sampler->sample_trilinear(tex,tx,ty,1.0f/(x1-x0),1.0f/(y1-y0)); 

          alpha_compost( &supersample_target[4 * ( (x + y * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] , color );
        }
}

void SoftwareRendererImp::rasterize_ellipse( float cx, float cy,
                                             float rx, float ry,
                                             Color color ) {
  // Task ?: 
  // Implement image rasterization

  for(int y = max((int)floor(cy-ry),0); y<=min((int)ceil(cy+ry),(int)target_h-1); ++y)
    for(int x = max((int)floor(cx-rx),0); x<=min((int)ceil(cx+rx),(int)target_w-1); ++x)
      for( int dy = 0; dy < sample_rate; ++dy )
        for( int dx = 0; dx < sample_rate; ++dx ) {
          float px = (float) x + 1.0 * (2*dx+1) / (2*sample_rate), py = (float) y + 1.0 * (2*dy+1) / (2*sample_rate);
          if( ((px-cx)*(px-cx))/(rx*rx) + ((py-cy)*(py-cy))/(ry*ry) <= 1.0 ) {
            alpha_compost( &supersample_target[4 * ( (x + y * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] , color );
          }
        }
}


// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 3: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 3".

  for( int y = 0; y < target_h; ++y )
    for( int x = 0; x < target_w; ++x ) {
      Color color=Color(0,0,0,0);
      for( int dy = 0; dy < sample_rate; ++dy )
        for( int dx = 0; dx < sample_rate; ++dx ) {
          Color tmpColor;
          uint8_to_float( &tmpColor.r , &supersample_target[4 * ( (x + y * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] );
          tmpColor.r/=(float)(sample_rate*sample_rate);
          tmpColor.g/=(float)(sample_rate*sample_rate);
          tmpColor.b/=(float)(sample_rate*sample_rate);
          tmpColor.a/=(float)(sample_rate*sample_rate);
          color+=tmpColor;
        }
      float_to_uint8( &render_target[4 * (x + y * target_w)] , &color.r );
    }

  return;
}


} // namespace CMU462
