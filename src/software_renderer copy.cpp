#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <queue>
#include <utility>

#include "triangulation.h"

using namespace std;

namespace CMU462 {

const float eps=1.0e-20;
int xinc[4]={0,0,-1,1};
int yinc[4]={-1,1,0,0};

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

  //cerr<<svg.width<<" "<<svg.height<<"\n";

  supersample_target = new unsigned char[4 * sample_rate * sample_rate * target_w * target_h];
  memset(supersample_target, 255, 4 * sample_rate * sample_rate * target_w * target_h);
  //------ Added --------------
  targeted = new int[4 * sample_rate * sample_rate * target_w * target_h];
  memset(targeted, 0, 16 * sample_rate * sample_rate * target_w * target_h);

  dst = new int[4 * sample_rate * sample_rate * target_w * target_h];
  memset(dst, 0, 16 * sample_rate * sample_rate * target_w * target_h);
  //------ End Added ----------

  // -------- Removed -------------
  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    group_transformation = canvas_to_screen;
    draw_element(svg.elements[i]);
  }
  // -------- End Removed ---------

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

  // -------- Added -------------

  /*Vector2D A=transform(Vector2D(svg.width/2,5));
  Vector2D B=transform(Vector2D(5,svg.height-5));
  Vector2D C=transform(Vector2D(svg.width-5,svg.height-5));

  //rasterize_line(p1.x, p1.y, p2.x, p2.y, Color::Black);
  //rasterize_point(p1.x, p1.y, Color::Black);
  //rasterize_point(p1.x+1, p1.y, Color::Black);
  //rasterize_point(p1.x, p1.y+1, Color::Black);
  //rasterize_point(p1.x+1, p1.y+1, Color::Black);

  rasterize_line(A.x, A.y, B.x, B.y, Color(1.0,0.0,0.0));
  rasterize_line(B.x, B.y, C.x, C.y, Color(1.0,0.0,0.0));
  rasterize_line(C.x, C.y, A.x, A.y, Color(1.0,0.0,0.0));

  Vector2D p=(A+B+C)/3.0;
  //for(int i=0; i<1000000; ++i)
  for(int i=0; i<1000000; ++i)
  {
    int r=std::rand()%3;
    if(r==0) p=(p+A)/2;
    else if(r==1) p=(p+B)/2;
    else p=(p+C)/2;
    rasterize_point(p.x, p.y, Color::Black);
    //rasterize_point(p.x+1, p.y, Color::Black);
    //rasterize_point(p.x, p.y+1, Color::Black);
    //rasterize_point(p.x+1, p.y+1, Color::Black);
  }*/

  // -------- End Added -------------

  // -------- Added -----------------

  /*Vector2D p(0.0,0.0);
  Vector2D tp=transform(Vector2D((p.x+3.5)*svg.width/7.0, (p.y+1)*svg.height/12.0));
  rasterize_point_2(tp.x, tp.y, Color(0.0,1.0,0.0));

  for(int i=0; i<1000000; ++i) {
    Vector2D q;

    int r=std::rand()%100;
    if(r<1) {
      q.x=0;
      q.y=0.16*p.y;
    }
    else if(r<86) {
      q.x=0.85*p.x+0.04*p.y;
      q.y=-0.04*p.x+0.85*p.y+1.6;
    }
    else if(r<93) {
      q.x=0.2*p.x-0.26*p.y;
      q.y=0.23*p.x+0.22*p.y+1.6;
    }
    else {
      q.x=-0.15*p.x+0.28*p.y;
      q.y=0.26*p.x+0.24*p.y+0.44;
    }

    p=q;
    Vector2D tp=transform(Vector2D((p.x+3.5)*svg.width/7.0, (p.y+1.0)*svg.height/12.0));
    rasterize_point_2(tp.x, tp.y, Color(0.0,1.0,0.0));
    //rasterize_point(tp.x, tp.y, Color::Black);
  }*/

  // -------- End Added -------------  

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
    case IFS:
      draw_ifs(static_cast<Ifs&>(*element));
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

  //delete[] supersample_target;
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

void SoftwareRendererImp::draw_ifs( Ifs& ifs ) {

  //cerr<<ifs.transformations[0]<<ifs.probabilities[0]<<"\n";
  /*Vector2D p(0.0,0.0);
  Vector2D tp=transform(Vector2D((p.x+3.5)*600.0/7.0, (p.y+1.0)*600.0/12.0));
  rasterize_point_2(tp.x, tp.y, Color(0.0,1.0,0.0));

  for(int i=0; i<1000000; ++i) {
    Vector2D q;

    float r=(float)(std::rand()) / RAND_MAX;
    if(r<=0.01) {
      q.x=0;
      q.y=0.16*p.y;
    }
    else if(r<=0.86) {
      q.x=0.85*p.x+0.04*p.y;
      q.y=-0.04*p.x+0.85*p.y+1.6;
    }
    else if(r<=0.93) {
      q.x=0.2*p.x-0.26*p.y;
      q.y=0.23*p.x+0.22*p.y+1.6;
    }
    else {
      q.x=-0.15*p.x+0.28*p.y;
      q.y=0.26*p.x+0.24*p.y+0.44;
    }

    p=q;
    if(i<10) cout<<p<<"\n";
    Vector2D tp=transform(Vector2D((p.x+3.5)*600.0/7.0, (p.y+1.0)*600.0/12.0));
    rasterize_point_2(tp.x, tp.y, Color(0.0,1.0,0.0));
    //rasterize_point(tp.x, tp.y, Color::Black);
  }*/

  /*for(int i=0; i<4; ++i)
    cerr<<ifs.transformations[i]<<ifs.probabilities[i]<<"\n";*/



    //adding sepenski
    if(!ifs.type){
      Vector2D c(300,300);
      Vector2D cn = transform(c);
      rasterize_point_2(cn.x,cn.y,Color::Black);
      /*for(int i=0;i<100000;i++)
      {
        int randNum = rand()%ifs.num_points;
        Vector2D dir = -transform(ifs.points[randNum]) + c;
         c = (transform(ifs.points[randNum]) + dir)/(ifs.r);
         //cout << randNum << " : " <<c << endl;
          rasterize_point_2(c.x,c.y,Color::Black);
      }*/
      /*ifs.points[0] = transform(ifs.points[0]);
      ifs.points[1] = transform(ifs.points[1]);
      ifs.points[2] = transform(ifs.points[2]);*/
      /*rasterize_line(ifs.points[0].x, ifs.points[0].y, ifs.points[1].x, ifs.points[1].y, Color::Black);
      rasterize_line(ifs.points[1].x, ifs.points[1].y, ifs.points[2].x, ifs.points[2].y, Color::Black);
      rasterize_line(ifs.points[2].x, ifs.points[2].y, ifs.points[0].x, ifs.points[0].y, Color::Black);
      rasterize_point_2(ifs.points[0].x,ifs.points[0].y,Color(1.0,0.0,0.0));
      rasterize_point_2(ifs.points[1].x,ifs.points[1].y,Color(1.0,0.0,0.0));
      rasterize_point_2(ifs.points[2].x,ifs.points[2].y,Color(1.0,0.0,0.0));

      cout << "p0: " << ifs.points[0] << endl;
      cout << "p1: " << ifs.points[1] << endl;
      cout << "p2: " << ifs.points[2] << endl;*/
      for(int i=0;i<1000000;i++)
      {
        int randNum = rand()%ifs.num_points;
        //double r = 1/(2*(1+cos(2*3.14159265/5)));
        //cout << r << endl;
        Vector2D dir = -ifs.points[randNum] + c;
         c.x = ifs.points[randNum].x + dir.x*(ifs.r);
         c.y = ifs.points[randNum].y + dir.y*(ifs.r);
         //cout << randNum << " : " <<c << endl;
         Vector2D cn = transform(c);
        rasterize_point_2(cn.x,cn.y,Color::Black);
          //rasterize_point(c.x,c.y,Color::Black);
      }
      return;

    }


  Vector3D p(ifs.seed.x,ifs.seed.y,1.0);
  Vector3D tp=ifs.renderTransformation*p;
  Vector2D rp=transform(Vector2D(tp.x,tp.y));
  rasterize_point_2(rp.x, rp.y, Color(0.0,1.0,0.0));

  for(int iter=0; iter<1000000; ++iter) {
    float r=(float)(std::rand()) / RAND_MAX;

    float cp=0.0;
    for(int i=0; i<(int)ifs.transformations.size(); ++i) {
      cp+=ifs.probabilities[i];
      //cerr<<cp<<"\n";

      if(r<=cp) {
        //if(iter<10) cout<<"@ "<<i<"\n";
        p = ifs.transformations[i] * p;
        p.x = sin(p.x)*cos(p.y);
        p.y = tan(p.y);


        /*
        if(i==1)
        {
          p.x = sin(p.x);
          p.y = sin(p.y);
        }
        else if(i==2)
        {
          double radiusSquared = (p.x*p.x+p.y*p.y);
          p.x = p.x/radiusSquared;
          p.y = p.y/radiusSquared;
        }
        else
        {
          double radiusSquared = (p.x*p.x+p.y*p.y);
          double x_old = p.x;
          double y_old = p.y;
          p.x = x_old*sin(radiusSquared) - y_old*cos(radiusSquared);
          p.y = x_old*cos(radiusSquared) + y_old*sin(radiusSquared);
        }*/
        //p = 
        break;
      }
    }
    //if(iter<10) cout<<"\n"<<p<<"\n";

    tp=ifs.renderTransformation*p;
    tp.x = 1.2*tp.x + 0;
    tp.y = 1.2*tp.y + 100;
    Vector2D rp=transform(Vector2D(tp.x,tp.y));
    rasterize_point_2(rp.x, rp.y, Color(0.0,1.0,0.0));
  }
  //cout<<ifs.renderTransformation;
}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

// -------- Added ----------------
void SoftwareRendererImp::add_to_target(int ix, int iy, Color color) {
  int sx=ix/sample_rate, dx=ix%sample_rate;
  int sy=iy/sample_rate, dy=iy%sample_rate;

  unsigned char tmp[4];
  float_to_uint8( tmp , &color.r );

  /*if(color.r>0.0)
    cerr<<" @@@ "<<color
      <<" --- "<<(int)supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx ) + 0]
      <<" --- "<<(int)supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx ) + 1]
      <<" --- "<<(int)supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx ) + 2]<<"\n"
      <<" --- "<<(int)tmp[0]
      <<" --- "<<(int)tmp[1]
      <<" --- "<<(int)tmp[2]<<"\n";*/

  if(targeted[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] == 0)
  {
      targeted[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] = 1;
      supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] = 0;
      supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx) + 1] = 0;
      supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx) + 2] = 0;
  }

  supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] =
    (unsigned char)min((int)(supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )])+(int)tmp[0],255);
  supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx) + 1] =
    (unsigned char)min((int)(supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx) + 1])+(int)tmp[1],255);
  supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx) + 2] =
    (unsigned char)min((int)(supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx) + 2])+(int)tmp[2],255);
  supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx ) + 3] = 255;
  //supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx ) + 3] =
    //min(supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx ) + 3] , (unsigned char)255);

  rmax=max(rmax,tmp[0]/255.0f);
  gmax=max(gmax,tmp[1]/255.0f);
  bmax=max(bmax,tmp[2]/255.0f);
}
// -------- End Added ------------

// -------- Added ---------------
void SoftwareRendererImp::rasterize_point_2( float x, float y, Color color ) {

  if(x<0 or x>target_w) return;
  if(y<0 or y>target_h) return;

  x = max( x , 1.0f/(2.0f*sample_rate) + eps ); x = min( x , target_w - 1.0f/(2.0f*sample_rate) - eps );
  y = max( y , 1.0f/(2.0f*sample_rate) + eps ); y = min( y , target_h - 1.0f/(2.0f*sample_rate) - eps );

  float ex=x*sample_rate;
  float ey=y*sample_rate;

  int ix=(int)(ex-0.5);
  int iy=(int)(ey-0.5);

  float fx=ix+0.5;
  float fy=iy+0.5;

  float x1=ex-fx, y1=ey-fy;
  float x2=1.0-x1, y2=1.0-y1;

  add_to_target(ix,iy,x2*y2*color);
  add_to_target(ix,iy+1,x2*y1*color);
  add_to_target(ix+1,iy,x1*y2*color);
  add_to_target(ix+1,iy+1,x1*y1*color);
}
// -------- End Added -----------

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

  // -------- Added ---------------

  /*float rmax=0.0,gmax=0.0,bmax=0.0;
  for( int y = 0; y < target_h; ++y ) {
    for( int x = 0; x < target_w; ++x ) {  
      for( int dy = 0; dy < sample_rate; ++dy ) {
        for( int dx = 0; dx < sample_rate; ++dx ) {
          Color tmpColor;
          uint8_to_float( &tmpColor.r , &supersample_target[4 * ( (x + y * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] );
          rmax=max(rmax,tmpColor.r);
          gmax=max(gmax,tmpColor.g);
          bmax=max(bmax,tmpColor.b);
        }
      }
    }
  }*/

  rmax=max(rmax,eps);
  gmax=max(gmax,eps);
  bmax=max(bmax,eps);

  for( int y = 0; y < target_h; ++y ) {
    for( int x = 0; x < target_w; ++x ) {
      for( int dy = 0; dy < sample_rate; ++dy ) {
        for( int dx = 0; dx < sample_rate; ++dx ) {
          Color tmpColor;
          uint8_to_float( &tmpColor.r , &supersample_target[4 * ( (x + y * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] );
          
          tmpColor.r/=rmax;
          tmpColor.g/=gmax;
          tmpColor.b/=bmax;

          float_to_uint8( &supersample_target[4 * ( (x + y * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] , &tmpColor.r );
        }
      }
    }
  }

  //cerr<<rmax<<" "<<gmax<<" "<<bmax<<"\n";

  queue<pair<int,int> > Q;

  for(int ix=0; ix<sample_rate*target_w; ++ix)
    for(int iy=0; iy<sample_rate*target_h; ++iy) {
      int sx=ix/sample_rate, dx=ix%sample_rate;
      int sy=iy/sample_rate, dy=iy%sample_rate;
      if(targeted[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )]==1)
      {
        Q.push(pair<int,int>(ix,iy));
        dst[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )]=0;
      }
    }

  int maxdst=0;
  while(not Q.empty()) {
    pair<int,int> cpt=Q.front(); Q.pop();
    int ix=cpt.first, iy=cpt.second;
    int sx=ix/sample_rate, dx=ix%sample_rate;
    int sy=iy/sample_rate, dy=iy%sample_rate;

    for(int dir=0; dir<4; ++dir)
    {
      int nix=ix+xinc[dir], niy=iy+yinc[dir];
      if(0<=nix and nix<sample_rate*target_w and 0<=niy and niy<sample_rate*target_h) {
        int nsx=nix/sample_rate, ndx=nix%sample_rate;
        int nsy=niy/sample_rate, ndy=niy%sample_rate;

        if(not targeted[4 * ( (nsx + nsy * target_w) * sample_rate * sample_rate + sample_rate * ndy + ndx )]) {
          dst[4 * ( (nsx + nsy * target_w) * sample_rate * sample_rate + sample_rate * ndy + ndx )] =
            dst[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] + 1;
          maxdst=max(maxdst,dst[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )] + 1);
          targeted[4 * ( (nsx + nsy * target_w) * sample_rate * sample_rate + sample_rate * ndy + ndx )] = 1;
          Q.push(pair<int,int>(nix,niy));
        }
      }
    }
  }

  for( int sy = 0; sy < target_h; ++sy ) {
    for( int sx = 0; sx < target_w; ++sx ) {
      for( int dy = 0; dy < sample_rate; ++dy ) {
        for( int dx = 0; dx < sample_rate; ++dx ) {
          if(dst[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )]>0)
          {
            Color tmpColor=Color(0.0,1.0,0.0)*(1.0*(dst[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )]-1)/(maxdst-1));
            float_to_uint8(&supersample_target[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )],&tmpColor.r);

            //cerr<<sx<<" "<<sy<<" : "<<dst[4 * ( (sx + sy * target_w) * sample_rate * sample_rate + sample_rate * dy + dx )]<<" "<<maxdst<<" "<<tmpColor<<"\n";
          }
        }
      }
    }
  }

  // -------- End Added -----------

  for( int y = 0; y < target_h; ++y ) {
    for( int x = 0; x < target_w; ++x ) {
      Color color=Color(0,0,0,0);
      for( int dy = 0; dy < sample_rate; ++dy ) {
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
    }
  }

  return;
}


} // namespace CMU462
