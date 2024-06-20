#ifndef RD_ENGINE_DIRECT_H
#define RD_ENGINE_DIRECT_H

#include "rd_enginebase.h"

#include <string>
using std::string;

// This is a rendering engine that renders directly to the image buffer as
// primitives come in.  A depth buffer is obviously needed.  Transparency is
// not properly handled.

class REDirect: public RenderEngine
{
 public:
  // Only methods inherited from the RenderEngine class should be added here,
  // as needed.
  virtual int rd_display(const string & name, const string & type,
			 const string & mode) override;

  virtual int rd_format(int xresolution, int yresolution) override;

  virtual int rd_world_begin(void) override;
  virtual int rd_world_end(void) override;

  virtual int rd_frame_begin(int frame_no) override;
  virtual int rd_frame_end(void) override;

  virtual int rd_render_cleanup(void) override;


  virtual int rd_circle(const float center[3], float radius) override;

  virtual int rd_line(const float start[3], const float end[3]) override;

  virtual int rd_point(const float p[3]) override;

   virtual int rd_background(const float color[]) override;
  // red, green, blue by default

  virtual int rd_color(const float color[]) override;

  virtual int rd_fill(const float seed_point[3]) override;


  virtual int rd_camera_eye(const float eyepoint[3]) override;
  virtual int rd_camera_at(const float atpoint[3]) override;
  virtual int rd_camera_up(const float up[3]) override;
  virtual int rd_camera_fov(float fov) override;
  virtual int rd_clipping(float znear, float zfar) override;

  virtual int rd_translate(const float offset[3]) override;
  virtual int rd_scale(const float scale_factor[3]) override;
  virtual int rd_rotate_xy(float angle) override;
  virtual int rd_rotate_yz(float angle) override;
  virtual int rd_rotate_zx(float angle) override;

  virtual int rd_pointset(const string & vertex_type,
			  int nvertex, const vector<float> & vertex) override;
  virtual int rd_polyset(const string & vertex_type,
			 int nvertex, const vector<float> & vertex,
			 int nface,   const vector<int> & face) override;

  virtual int rd_cone(float height, float radius, float thetamax) override;
  virtual int rd_cube(void) override;
  virtual int rd_cylinder(float radius, float zmin,
			  float zmax, float thetamax) override;
  virtual int rd_disk(float height, float radius, float theta) override;

  virtual int rd_sphere(float radius, float zmin, float zmax, float thetamax) override;


  virtual int rd_xform_push(void) override;
  virtual int rd_xform_pop(void) override;

};

#endif /* RD_ENGINE_DIRECT_H */
