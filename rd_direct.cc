#include <iostream>
#include <math.h>
#include<queue>
#include<vector>
#include "rd_direct.h"
#include "rd_enginebase.h"
#include "pnm_display.h"
#include "rd_error.h"
#include "rd_display.h"
#include "pp_transform.h"
#include "attr_point.h"

#include <string>
using std::string;
using namespace std;

//global vriables
int frame_number = 0;

float drawing_color [] = {1.0,1.0,1.0};

point camera_eye = {0.0,0.0,0.0};
point camera_at = {0.0, 0.0,-1.0};
Vector camera_up = {0.0,0.1,0.0};
float camera_fov = 90.0;
float near = 1.0, far = 1.0e+09;
pointh old_vertex;

float** depth_buffer = nullptr;

Edge** edge_table = nullptr;
bool edge_table_initialised = false;

xform current_transform;
xform world_to_clipping_matrix;
xform clipping_to_divice_matrix;

void rd_draw_line(int x0, int x1, int y0, int y1);
void rd_draw_circle(int x, int y, int radius);
void rd_flood_fill(int x, int y, float currentColor[3]);
void circle_using_angle(float x, float y, float z,const float radius);
void line_pipeline(const pointh& vertex, string flag);
void point_pipeline(const pointh& vertex);
void clip_line(pointh& start, pointh& end);
void cube();
void rd_display_pixel(int x, int y, float z);
int poly_pipeline(attr_point p, int end_flag);
void scan_convert(attr_point* clipped_list, int count);
bool buildEdgeList(const attr_point* points, int count);
void addActiveList(int scanline, Edge*& AET);
void makeEdgeRec(attr_point upper, attr_point lower);
void insertEdge(Edge* &list, Edge* e);
void updateAET(int scanline, Edge*& active);
void deleteAfter(Edge* q);
void resortAET(Edge*& active);
void fillBetweenEdges(int scanline, Edge* AET);
int poly_clip(int n_vertex, attr_point vertex_list[], attr_point clipped_list[]);
int calculate_clip_code(const pointh& point);
bool inside(const attr_point& p , Boundary b);
void clip_a_point(const attr_point& p, Boundary b, attr_point (&first)[], attr_point (&last)[], bool (&flag)[], attr_point clipped_list[], int& count);
void clip_last_point(attr_point (&first)[],attr_point (&last)[], bool (&flag)[], attr_point clipped_list[], int& count);
float calculate_intersection(const pointh& start, const pointh& end, int clip_code) ;
attr_point intersect(const attr_point& p1, const attr_point& p2, Boundary b);
pointh interpolate_point(const pointh& start, const pointh& end, float t) ;
bool cross(const attr_point& p1, const attr_point& p2, Boundary b);
//implementing RenderEngine methods
int REDirect::rd_display(const string & name, const string & type,
			 const string & mode){

    return RD_OK;
}

int REDirect::rd_format(int xresolution, int yresolution){

    display_xSize = xresolution;
    display_ySize = yresolution;
    return RD_OK;
}

int REDirect::rd_camera_eye(const float eyepoint[3]){

    camera_eye.x=eyepoint[0];
    camera_eye.y=eyepoint[1];
    camera_eye.z=eyepoint[2];

    return RD_OK;
}

int REDirect::rd_camera_at(const float atpoint[3]){

    camera_at.x=atpoint[0];
    camera_at.y=atpoint[1];
    camera_at.z=atpoint[2];

    return RD_OK;
}

int REDirect::rd_camera_up(const float up[3]){

    camera_up.x=up[0];
    camera_up.y=up[1];
    camera_up.z=up[2];

    return RD_OK;
}

int REDirect::rd_camera_fov(float fov){

    camera_fov = fov;

    return RD_OK;
}

int REDirect::rd_clipping(float znear, float zfar){

    near = znear;
    far = zfar;

    return RD_OK;
}

int REDirect::rd_frame_begin(int frame_no){

    frame_number = frame_no;
    return RD_OK;
}

int REDirect::rd_world_begin(){



    if(!depth_buffer){

        depth_buffer = new float*[display_ySize];
        for(int i=0;i<display_ySize;i++){
            depth_buffer[i] = new float[display_xSize];

            for (int j = 0; j < display_xSize; ++j) {
                depth_buffer[i][j] = std::numeric_limits<float>::max();
            }
        }
    }


        edge_table = new Edge*[display_ySize];
        for (int i=0;i<display_ySize;i++){
            edge_table[i] = new Edge();


        }



    rd_disp_init_frame(frame_number);
    current_transform = identity();

    xform world_to_camera_matrix = world_to_camera(camera_eye,camera_at,camera_up);

    double aspect_ratio = static_cast<double>(display_xSize)/static_cast<double>(display_ySize);
    xform camera_to_clipping_matrix = camera_to_clip(camera_fov, near, far,aspect_ratio);

    world_to_clipping_matrix = multiply(camera_to_clipping_matrix,world_to_camera_matrix);

    clipping_to_divice_matrix = clip_to_device(display_xSize, display_ySize);

    return RD_OK;
}

int REDirect::rd_world_end(){

    if(depth_buffer){
        for(int i = 0; i < display_ySize; i++){
            for(int j = 0; j < display_xSize; j++){
                depth_buffer[i][j] = std::numeric_limits<float>::max();
            }
        }
    }

    if(edge_table){
        for (int i = 0; i < display_ySize; ++i) {
            delete edge_table[i];
        }
        delete[] edge_table;
    }

    if (rd_disp_end_frame() != RD_OK) {
        return RD_INPUT_EXPECTED_EOF;
    }
    return RD_OK;
}

int REDirect::rd_render_cleanup(){

    if(depth_buffer){
        for(int i=0;i<display_ySize;i++){
            delete[] depth_buffer[i];
        }
        delete[] depth_buffer;
        depth_buffer=nullptr;
    }

    return RD_OK;
}

int REDirect::rd_frame_end(){
    if (rd_disp_end_frame() != RD_OK) {
        return RD_INPUT_EXPECTED_EOF;
    }
    return RD_OK;
}

int REDirect::rd_translate(const float offset[3]){

    float tx = offset[0];
    float ty = offset[1];
    float tz = offset[2];

    current_transform = translate(current_transform, tx, ty, tz);
    return RD_OK;
}

int REDirect::rd_scale(const float scale_factor[3]){

    float sx = scale_factor[0];
    float sy = scale_factor[1];
    float sz = scale_factor[2];

    current_transform = scale(current_transform, sx, sy, sz);
    return RD_OK;
}

int REDirect::rd_rotate_xy(float angle){

    rotate_xy(current_transform, angle);
    return RD_OK;
}

int REDirect::rd_rotate_yz(float angle){

    rotate_yz(current_transform, angle);
    return RD_OK;
}

int REDirect::rd_rotate_zx(float angle){

    rotate_zx(current_transform, angle);
    return RD_OK;
}

int REDirect::rd_xform_push(){

    push(current_transform);
    return RD_OK;
}

int REDirect::rd_xform_pop(){

    current_transform = pop();
    return RD_OK;
}
int REDirect::rd_color(const float color[]){

    if(color!=nullptr){
        for(int i=0;i<3;i++){
            drawing_color[i] = color[i];
        }
    }

    return RD_OK;
}

int REDirect::rd_background(const float color[]){


    if(color!=nullptr){
        rd_set_background(color);
    }
    else{
        const float black_color[] = {0.0,0.0,0.0};
        rd_set_background(black_color);
    }

    return RD_OK;
}

int REDirect::rd_point(const float p[3]){

   /* if(p==nullptr){
        return RD_INPUT_UNINITIALIZED_DISPLAY;
    }

    int x = static_cast<int>(p[0]);
    int y = static_cast<int>(p[1]);

    if(x<0 || y<0 || x>display_xSize || y>display_ySize){
        return RD_OK;
    }
    else{
        rd_write_pixel(x, y, drawing_color);
    }
    */
    point start_point;
    start_point.x = p[0];
    start_point.y = p[1];
    start_point.z = p[2];

    pointh start_pointh = convert(start_point);
    point_pipeline(start_pointh);
    return RD_OK;
}

int REDirect::rd_line(const float start[3], const float end[3]){
/*
    int x0 = static_cast<int>(start[0]);
    int x1 = static_cast<int>(end[0]);

    int y0 = static_cast<int>(start[1]);
    int y1 = static_cast<int>(end[1]);



    rd_draw_line(x0,x1,y0,y1);
    */

    point start_point;
    start_point.x = start[0];
    start_point.y = start[1];
    start_point.z = start[2];

    point end_point;
    end_point.x = end[0];
    end_point.y = end[1];
    end_point.z = end[2];

    pointh start_pointh = convert(start_point);
    pointh end_pointh = convert(end_point);

    line_pipeline(start_pointh, "MOVE");
    line_pipeline(end_pointh, "DRAW");
    return RD_OK;
}

int REDirect::rd_circle(const float center[3], float radius){

    int x = static_cast<int>(center[0]);
    int y = static_cast<int>(center[1]);

    rd_draw_circle(x, y, static_cast<int>(radius));

    return RD_OK;
}

int REDirect::rd_fill(const float seed_point[3]){

    int x = static_cast<int>(seed_point[0]);
    int y = static_cast<int>(seed_point[1]);

    if(!(x>=0 && x < display_xSize && y>=0 && y< display_ySize))
    {
        return RD_INPUT_ILLEGAL_VERTEX_INDEX;
    }

    float currentColor[3];
    rd_read_pixel(x, y, currentColor);

    rd_flood_fill(x, y, currentColor);

    return RD_OK;
}

int REDirect::rd_pointset(const string & vertex_type, int nvertex, const vector<float> & vertex){

    for(size_t i=0;i<vertex.size();i=i+3){
        point p;
        p.x = vertex[i];
        p.y=vertex[i+1];
        p.z = vertex[i+2];

        pointh ph = convert(p);

        point_pipeline(ph);
    }

    return RD_OK;
}

int REDirect::rd_polyset(const string & vertex_type, int nvertex, const vector<float> & vertex, int nface,   const vector<int> & face){


     vector<attr_point> vertices;
    for (size_t i = 0; i < vertex.size(); i += 3) {
        attr_point p = { vertex[i], vertex[i + 1], vertex[i + 2], 1.0f };
        vertices.push_back(p);
    }

    vector<int> face_array;

    for (size_t i = 0; i < face.size(); i++) {
        if (face[i] == -1) {
            for (size_t k = 0; k < face_array.size(); k++) {
                poly_pipeline(vertices[face_array[k]], (k == face_array.size() - 1) ? 1 : 0);
            }
            face_array.clear();
        } else {
            face_array.push_back(face[i]);
        }
    }

    return RD_OK;
}
/*
int REDirect::rd_cone(float height, float radius, float thetamax){

    float x = radius,y=0,z=0;


    circle_using_angle(x,y,z,radius);

    const int NSTEPS = 20;
    for (int i = 0; i <= NSTEPS; ++i) {

        float theta = i * (2 * M_PI) / NSTEPS;


        float x = radius * cos(theta);
        float y = radius * sin(theta);

        pointh bottom_point = { x, y, 0.0f, 1.0f };
        pointh top_point = { 0.0f,0.0f, height, 1.0f };
        line_pipeline(bottom_point, "MOVE");
        line_pipeline(top_point, "DRAW");
    }
    return RD_OK;
}
*/
int REDirect::rd_cone(float height, float radius, float thetamax) {
    const int NSTEPS = 20;
    const float step_theta = (2 * M_PI) / NSTEPS;


    attr_point top = { 0.0f, 0.0f, height, 1.0f };

    for (int i = 0; i < NSTEPS; ++i) {
        float theta = i * step_theta;
        float theta1 = (i + 1) * step_theta;

        float x1 = radius * cos(theta);
        float y1 = radius * sin(theta);
        float x2 = radius * cos(theta1);
        float y2 = radius * sin(theta1);

        attr_point bottom1 = { x1, y1, 0.0f, 1.0f };
        attr_point bottom2 = { x2, y2, 0.0f, 1.0f };

        vector<attr_point> vertices;
        vertices.push_back(top);
        vertices.push_back(bottom1);
        vertices.push_back(bottom2);

        for(int j=0;j<3;j++){
            poly_pipeline(vertices[j], (j==2)?1:0);
        }

    }

    return RD_OK;
}
int REDirect::rd_cube(){
   // cube();
   const int cube_vertex_count = 8;
    attr_point cube_vertices[cube_vertex_count] = {
        {{-1.0f, -1.0f, -1.0f, 1.0f}}, // Bottom front left
        {{1.0f, -1.0f, -1.0f, 1.0f}},  // Bottom front right
        {{1.0f, 1.0f, -1.0f, 1.0f}},   // Top front right
        {{-1.0f, 1.0f, -1.0f, 1.0f}},  // Top front left
        {{-1.0f, -1.0f, 1.0f, 1.0f}},  // Bottom back left
        {{1.0f, -1.0f, 1.0f, 1.0f}},   // Bottom back right
        {{1.0f, 1.0f, 1.0f, 1.0f}},    // Top back right
        {{-1.0f, 1.0f, 1.0f, 1.0f}}    // Top back left
    };

    const int cube_face_indices[][4] = {
       {2, 3, 0, 1},  // Front face
    {6, 2, 1, 5},  // Right face
        {4, 5, 6, 7},  // Back face
        {0, 4, 7, 3},  // Left face
       {0, 1, 5, 4},  // Bottom face
        {3, 2, 6, 7}   // Top face
    };

    const int num_faces = sizeof(cube_face_indices) / sizeof(cube_face_indices[0]);
    for (int i = 0; i < num_faces; ++i) {
       const int* face_indices = cube_face_indices[i];
        for (int j = 0; j < 4; ++j) {
            poly_pipeline(cube_vertices[face_indices[j]], (j == 3) ? 1 : 0);
        }
    }


    return RD_OK;
}
/*
int REDirect::rd_cylinder(float radius, float zmin, float zmax, float thetamax){

    rd_disk(zmin, radius, thetamax);
    rd_disk(zmax, radius, thetamax);

    const int NSTEPS = 30;
    for (int i = 0; i <= NSTEPS; ++i) {

        float theta = i * (2 * M_PI) / NSTEPS;


        float x = radius * cos(theta);
        float y = radius * sin(theta);

        pointh bottom_point = { x, y, zmin, 1.0f };
        pointh top_point = { x, y, zmax, 1.0f };
        line_pipeline(bottom_point, "MOVE");
        line_pipeline(top_point, "DRAW");
    }

    return RD_OK;
}
*/
int REDirect::rd_cylinder(float radius, float zmin, float zmax, float thetamax) {
    const int NSTEPS = 20;
    const float step_theta = (2 * M_PI) / NSTEPS;

    for (int i = 0; i < NSTEPS; i++) {
        float theta = i * step_theta;
        float theta1 = (i + 1) * step_theta;

        float x1 = radius * cos(theta);
        float y1 = radius * sin(theta);
        float x2 = radius * cos(theta1);
        float y2 = radius * sin(theta1);

        attr_point bottom1 = { x1, y1, zmin, 1.0f };
        attr_point bottom2 = { x2, y2, zmin, 1.0f };
        attr_point top1 = { x1, y1, zmax, 1.0f };
        attr_point top2 = { x2, y2, zmax, 1.0f };

        vector<attr_point> bottom_face = { bottom1, bottom2 };
        vector<attr_point> top_face = { top1, top2 };
    /*
        // Draw the bottom face
        for (int j = 0; j < 2; ++j) {
            poly_pipeline(bottom_face[j], (j == 1) ? 1 : 0);  // Move for the last vertex
        }

        // Draw the top face
        for (int j = 0; j < 2; ++j) {
            poly_pipeline(top_face[j], (j == 1) ? 1 : 0);  // Move for the last vertex
        }
*/
        // Draw the sides
        poly_pipeline(bottom1, 0);  // Move to the first vertex of the bottom face
        poly_pipeline(top1, 0);     // Draw to the first vertex of the top face
        poly_pipeline(top2, 0);     // Move to the second vertex of the top face
        poly_pipeline(bottom2, 1);  // Draw to the second vertex of the bottom face
    }
    return RD_OK;
}

int REDirect::rd_disk(float height, float radius, float theta){

    circle_using_angle(0.0f, 0.0f, height, radius);
    return RD_OK;
}

int REDirect::rd_sphere(float radius, float zmin, float zmax, float thetamax){

    int NSTEPS = 15;

    //latitude (theta)
    for(int i = 0; i < NSTEPS; i++){

            float theta1 = static_cast<float>(i) / NSTEPS * M_PI;
            float theta2 = static_cast<float>(i + 1) / NSTEPS * M_PI;

        //longitude (phi)
        for(int j = 0; j < NSTEPS; j++){


            float phi1 = static_cast<float>(j) / NSTEPS * 2 * M_PI;
            float phi2 = static_cast<float>(j + 1) / NSTEPS * 2 * M_PI;

            // Calculate vertices
            float x1 = radius * sin(theta1) * cos(phi1);
            float y1 = radius * sin(theta1) * sin(phi1);
            float z1 = radius * cos(theta1);

            float x2 = radius * sin(theta1) * cos(phi2);
            float y2 = radius * sin(theta1) * sin(phi2);
            float z2 = radius * cos(theta1);

            float x3 = radius * sin(theta2) * cos(phi2);
            float y3 = radius * sin(theta2) * sin(phi2);
            float z3 = radius * cos(theta2);

            float x4 = radius * sin(theta2) * cos(phi1);
            float y4 = radius * sin(theta2) * sin(phi1);
            float z4 = radius * cos(theta2);


            attr_point p1 = {x1, y1, z1, 1.0f};
            attr_point p2 = {x2, y2, z2, 1.0f};
            attr_point p3 = {x3, y3, z3, 1.0f};
            attr_point p4 = {x4, y4, z4, 1.0f};


            //vector<attr_point> polygon = {vertex1, vertex2, vertex3, vertex4};
            poly_pipeline(p1, 0);
            poly_pipeline(p2, 0);
            poly_pipeline(p3, 0);
            poly_pipeline(p4, 1);

        }
    }

    return RD_OK;
}

//Bresenham's line algorithm

/*
void rd_draw_line(int x0, int x1, int y0, int y1){

    int dx = abs(x1-x0);
    int dy = -abs(y1-y0);

    //checking positive and negative directions
    int x_direction = (x1>x0)?1:-1;
    int y_direction = (y1>y0)?1:-1;


    if (dx == 0) {
        // Vertical line

        for (int y = min(y0, y1); y <= max(y0, y1); ++y) {
            rd_write_pixel(x0, y, drawing_color);
        }
        return;
    }

    int error = dx+dy;

    while(true){
        rd_write_pixel(x0,y0,drawing_color);

        if(x0 == x1 && y0==y1)
            break;

        int error2 = 2*error;

        if(error2 >= dy){
            if( x0 == x1)
                break;
            error = error + dy;
            x0 = x0+x_direction;

        }
        if(error2 <= dx){
            if (y0 == y1)
                break;
            error = error + dx;
            y0 = y0 + y_direction;
        }
    }

}
*/

/*
void rd_draw_line(int x0, int x1, int y0, int y1){

    int dx = x1-x0;
    int dy = y1-y0;
    int m = dy/dx;

    if (dx == 0) {
        // Vertical line

        for (int y = min(y0, y1); y <= max(y0, y1); ++y) {
            rd_write_pixel(x0, y, drawing_color);
        }
        return;
    }
    if(m<1){
        int p = 2*dy-dx;
        int x = x0;
        int y=y0;

        if(dx<0){
            x=x1;
            y=y1;
            x1=x0;
        }
        rd_write_pixel(x,y,drawing_color);

        while(x<x1){
            if(p>=0){
                x=x+1;
                y=y+1;
                p = p + 2*dy - 2*dx;
            }
            else{
                x=x+1;
                p = 2*dy;
            }
            rd_write_pixel(x,y,drawing_color);
        }
    }
    else if(m>1){

        int p = 2*dx - dy;
        int x=x0;
        int y=y0;
        if(dy<0){
            x=x1;
            y=y1;
            y1=y0;

        }
        rd_write_pixel(x,y,drawing_color);
        while(y<y1){
            if(p>=0){
                x++;
                y++;
                p = p + 2*dx - 2*dy;
            }
            else{
                y++;
                p = p+2*dx;
            }
            rd_write_pixel(x,y,drawing_color);
        }
    }
    else if(m==1){
        int x=x0;
        int y=y0;
        rd_write_pixel(x,y,drawing_color);
        while(x<x1){
            x++;
            y++;
            rd_write_pixel(x,y,drawing_color);
        }
    }
}
    //checking positive and negative directions
*/

/*    int x=x0, y=y0;

    int p = 2*dy - dx;

    for(int i=0; i<=dx; i++){

        rd_write_pixel(x,y,drawing_color);

        x += 1;

        if(p>0){
            y+= 1;
            p = p + 2*(dy-dx);
        }
        else{
            p += 2*dy;
        }
    }

*/

//DDA algorithm
void rd_draw_line(float x0, float y0, float z0, float x1, float y1, float z1){

    int dx = static_cast<int>(x1) - static_cast<int>(x0);
    int dy = static_cast<int>(y1) - static_cast<int>(y0);

    int NSTEPS = max(abs(dx), abs(dy));

    rd_display_pixel(static_cast<int>(x0), static_cast<int>(y0), (z0));
    for(int i=1;i<= NSTEPS; i++){
       float t = static_cast<float>(i) / static_cast<float>(NSTEPS);

        int x = static_cast<int>(x0 + t * (x1 - x0));
        int y = static_cast<int>(y0 + t * (y1 - y0));
        float z = (z0 + t * (z1 - z0));

        rd_display_pixel(x,y,z);
    }
}

void rd_draw_circle(int x1, int y1, int radius){

    int x=0, y = radius;

    int p = 1-radius;

    rd_write_pixel(x1 + x, y1 + y, drawing_color);
    rd_write_pixel(x1 - x, y1 + y, drawing_color);
    rd_write_pixel(x1 + x, y1 - y, drawing_color);
    rd_write_pixel(x1 - x, y1 - y, drawing_color);
    rd_write_pixel(x1 + y, y1 + x, drawing_color);
    rd_write_pixel(x1 - y, y1 + x, drawing_color);
    rd_write_pixel(x1 + y, y1 - x, drawing_color);
    rd_write_pixel(x1 - y, y1 - x, drawing_color);

    while(y>=x){


        x++;
        if(p>0){
            y--;
            p = p + 2*x + 1 - 2*y;

        }
        else{
            p = p + 2*x+1;
        }

        rd_write_pixel(x1 + x, y1 + y, drawing_color);
        rd_write_pixel(x1 - x, y1 + y, drawing_color);
        rd_write_pixel(x1 + x, y1 - y, drawing_color);
        rd_write_pixel(x1 - x, y1 - y, drawing_color);
        rd_write_pixel(x1 + y, y1 + x, drawing_color);
        rd_write_pixel(x1 - y, y1 + x, drawing_color);
        rd_write_pixel(x1 + y, y1 - x, drawing_color);
        rd_write_pixel(x1 - y, y1 - x, drawing_color);

    }
}

void rd_flood_fill(int x, int y, float color[3]){

    queue<pair<int, int>> points;

    points.push({x,y});

    while(!points.empty()){
        auto[x0, y0] = points.front();
        points.pop();

        if(x0 >= 0 && x0 < display_xSize && y0 >= 0 && y0 < display_ySize){
            float currentColor[3];
            rd_read_pixel(x0,y0, currentColor);

            if(currentColor[0] == color[0] && currentColor[1] == color[1] && currentColor[2]== color[2]){

                rd_write_pixel(x0,y0,drawing_color);

                points.push({x0+1, y0});
                points.push({x0-1, y0});
                points.push({x0, y0+1});
                points.push({x0, y0-1});

            }
        }
    }

}
static bool is_previous_clipped = false;
void line_pipeline(const pointh& vertex, string flag){

    //static bool is_previous_clipped = false;
    pointh world = multiply(current_transform, vertex);

    pointh clipping_point = multiply(world_to_clipping_matrix, world);

    if(flag == "MOVE"){
        old_vertex = clipping_point;
    }

    else if(flag == "DRAW"){

        pointh unclipped = clipping_point;

       if (!is_previous_clipped){
        clip_line(old_vertex,clipping_point);
       }


        is_previous_clipped = (flag == "MOVE" && clipping_point.w == 0.0f);

        pointh deviceVertex = multiply(clipping_to_divice_matrix, old_vertex);

        pointh currentVertex = multiply(clipping_to_divice_matrix, clipping_point);

        old_vertex = unclipped;

        if(deviceVertex.w != 0 && currentVertex.w != 0)
            rd_draw_line((deviceVertex.x/deviceVertex.w), (deviceVertex.y/deviceVertex.w), (deviceVertex.z/deviceVertex.w),(currentVertex.x/currentVertex.w),(currentVertex.y/currentVertex.w), (currentVertex.z/currentVertex.w));

    }
}

void point_pipeline(const pointh& vertex){

    pointh world = multiply(current_transform, vertex);

    pointh clipping_point = multiply(world_to_clipping_matrix, world);

    if(clipping_point.x/clipping_point.w>=0 && clipping_point.x/clipping_point.w <=1 && clipping_point.y/clipping_point.w>=0 && clipping_point.y/clipping_point.w <=1 && clipping_point.z/clipping_point.w>=0 && clipping_point.z/clipping_point.w <=1){

        pointh deviceVertex = multiply(clipping_to_divice_matrix, clipping_point);
        rd_display_pixel(deviceVertex.x/deviceVertex.w, deviceVertex.y/deviceVertex.w, deviceVertex.z/deviceVertex.w);
    }


}

void cube(){
    pointh vertices[8] = {
        {-1.0f, -1.0f, -1.0f, 1.0f},
        {1.0f, -1.0f, -1.0f, 1.0f},
        {1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, -1.0f, 1.0f, 1.0f},
        {1.0f, -1.0f, 1.0f, 1.0f},
        {1.0f, 1.0f, 1.0f, 1.0f},
        {-1.0f, 1.0f, 1.0f, 1.0f}
        };

        int faces[6][4] = {
            //{0,1,2,3},
            //{1,2,6,5}
          // {0,4,7,3},
          // {5,6,7,4}
           // {0,1,5,4},
           // {2,3,7,6}
        };

        for(int i=0;i<6;i++){
            line_pipeline(vertices[faces[i][0]], "MOVE");
            line_pipeline(vertices[faces[i][1]], "DRAW");
            line_pipeline(vertices[faces[i][2]], "DRAW");
            line_pipeline(vertices[faces[i][3]], "DRAW");
            //line_pipeline(vertices[faces[i][0]], "DRAW");
        }
}

void circle_using_angle(float x, float y, float z, const float radius){

    float w=1;
    x = radius;
    const int NSTEPS = 30;
    pointh ph = {x,y,z,w};
    line_pipeline(ph, "MOVE");

    for(int i=0;i<=NSTEPS;i++){


        float theta = static_cast<float>(i*2*M_PI/NSTEPS);

        x=radius*cos(theta);
        y=radius*sin(theta);

        ph = {x,y,z,w};

        line_pipeline(ph, "DRAW");
    }
}

void clip_line(pointh& start, pointh& end) {

    pointh start_temp = start;
    pointh end_temp = end;


    int kode0 = calculate_clip_code(start);
    int kode1 = calculate_clip_code(end);

    int BC0[6] = {static_cast<int>(start.x), static_cast<int>(start.w-start.x), static_cast<int>(start.y), static_cast<int>(start.w-start.y), static_cast<int>(start.z), static_cast<int>(start.w-start.z)};
    int BC1[6] = {static_cast<int>(end.x), static_cast<int>(end.w-end.x),static_cast<int>(end.y),static_cast<int>(end.w-end.y),static_cast<int>(end.z),static_cast<int>(end.w-end.z)};
    if (kode0 & kode1) {

        start = end = {0, 0, 0, 1};
        return;
    }

    if (kode0 | kode1) {
       // float t_start = calculate_intersection(BC0, BC1, kode0);
       // float t_end = calculate_intersection(BC0, BC1, kode1);

        int kode = kode0 | kode1;
        int mask = 1;
        float alpha0 = 0.0f, alpha1 = 1.0f;

        for(int i=0;i<6;i++, mask<<=1){
            if((kode & mask) == 0){
                continue;
            }
            else{
                float alpha = BC0[i]/(BC0[i]-BC1[i]);

                if(kode0 & mask){
                    //outside to inside
                    alpha0 = max(alpha0,alpha);
                }
                else{
                    //inside to outside
                    alpha1 = min(alpha1, alpha);
                }
                if(alpha1<alpha0){
                    break;
                }
            }

        }
        start = interpolate_point(start_temp, end_temp, alpha0);
        end = interpolate_point(start_temp, end_temp, alpha1);
    }
}

int calculate_clip_code(const pointh& point) {
    int clip_code = 0;

    if (point.x < 0) clip_code |= 1; // Set bit 1
    if ((point.w - point.x) < 0) clip_code |= 1 << 1; // Set bit 2
    if (point.y < 0) clip_code |= 1 << 2; // Set bit 3
    if ((point.w - point.y) < 0) clip_code |= 1 << 3; // Set bit 4
    if (point.z < 0) clip_code |= 1 << 4; // Set bit 5
    if ((point.w - point.z) < 0) clip_code |= 1 << 5; // Set bit 6
    return clip_code;
}

float calculate_intersection(const pointh& start, const pointh& end, int clip_code) {
    float t = 0.0;

    if(clip_code & 1<<1){
        t = start.x/(start.x - end.x);
    }
    else if(clip_code & 1<<2){
        t = (start.w-start.x)/((start.w-start.x) - (end.w-end.x));
    }
    else if(clip_code & 1<<3){
        t = start.y/(start.y - end.y);
    }
    else if(clip_code & 1<<4){
        t = (start.w-start.y)/((start.w-start.y) - (end.w-end.y));
    }
    else if(clip_code & 1<<5){
        t = start.z/(start.z - end.z);
    }
    else if(clip_code & 1<<6){
        t = (start.w-start.z)/((start.w-start.z) - (end.w-end.z));

    }
    return t;
}
pointh interpolate_point(const pointh& start, const pointh& end, float t) {
    pointh interpolated_point;
    interpolated_point.x = start.x + t * (end.x - start.x);
    interpolated_point.y = start.y + t * (end.y - start.y);
    interpolated_point.z = start.z + t * (end.z - start.z);
    interpolated_point.w = start.w;

    return interpolated_point;
}
void rd_display_pixel(int x, int y, float z){

        if(x<0 || x>=display_xSize || y<0 || y>= display_ySize){
            return;
    }

    if(z < depth_buffer[y][x]){
        depth_buffer[y][x] = z;

        rd_write_pixel(x,y,drawing_color);
    }
}
const int MAX_VERTEX_LIST_SIZE = 50;
static attr_point clipped_list[MAX_VERTEX_LIST_SIZE];
int poly_pipeline(attr_point p, int end_flag){

    pointh geom, norm;
    pointh dev;

    //const int MAX_VERTEX_LIST_SIZE = 50;
    static attr_point vertex_list[MAX_VERTEX_LIST_SIZE];
    //static attr_point clipped_list[MAX_VERTEX_LIST_SIZE];
    static int n_vertex = 0;

    int i;

    geom.x = p.coord[0];
    geom.y = p.coord[1];
    geom.z = p.coord[2];
    geom.w = p.coord[3];

    geom = multiply(current_transform , geom);

    geom = multiply(world_to_clipping_matrix , geom);

    p.coord[0] = geom.x;
    p.coord[1] = geom.y;
    p.coord[2] = geom.z;
    p.coord[3] = geom.w;

    if(n_vertex == MAX_VERTEX_LIST_SIZE)
        return -1;

    vertex_list[n_vertex] = p;
    n_vertex++;

    if(end_flag ==0) //MOVE
        return 0;

    if((n_vertex = poly_clip(n_vertex, vertex_list, clipped_list))){

        for(i=0;i<n_vertex; i++){

            dev.x = clipped_list[i].coord[0];
            dev.y = clipped_list[i].coord[1];
            dev.z = clipped_list[i].coord[2];
            dev.w = clipped_list[i].coord[3];
/*
            norm.x = clipped_list[i+1].coord[0];
            norm.y = clipped_list[i+1].coord[1];
            norm.z = clipped_list[i+1].coord[2];
            norm.w = clipped_list[i+1].coord[3];

            dev = multiply(clipping_to_divice_matrix, dev);
            norm = multiply(clipping_to_divice_matrix, norm);

*/
            dev = multiply(clipping_to_divice_matrix, dev);

            clipped_list[i].coord[0] = dev.x;
            clipped_list[i].coord[1] = dev.y;
            clipped_list[i].coord[2] = dev.z;
            clipped_list[i].coord[3] = dev.w;


            clipped_list[i].coord[0] /= clipped_list[i].coord[3];
            clipped_list[i].coord[1] /= clipped_list[i].coord[3];
            clipped_list[i].coord[2] /= clipped_list[i].coord[3];


           // rd_draw_line(dev.x/dev.w, dev.y/dev.w, dev.z/dev.w, norm.x/norm.w, norm.y/norm.w, norm.z/norm.w);
           // rd_display_pixel(clipped_list[i].coord[0], clipped_list[i].coord[1], clipped_list[i].coord[2]);

        }
        for(i=0;i<n_vertex-1;i++){

            rd_draw_line(clipped_list[i].coord[0], clipped_list[i].coord[1], clipped_list[i].coord[2], clipped_list[i+1].coord[0], clipped_list[i+1].coord[1], clipped_list[i+1].coord[2]);
        }
        //rd_draw_line(clipped_list[i].coord[0], clipped_list[i].coord[1], clipped_list[i].coord[2], clipped_list[0].coord[0], clipped_list[0].coord[1], clipped_list[0].coord[2]);

        scan_convert(clipped_list, n_vertex);
    }

    n_vertex = 0;

    return RD_OK;
}
void scan_convert(attr_point* clipped_list, int count) {

    Edge* AET = nullptr; //Active Edge Table

    if(!buildEdgeList(clipped_list, count)){
        return; //no edges cross a scanline
    }

    for(int scan=0;scan < display_ySize; scan++){

        addActiveList(scan, AET);

        if(AET!= nullptr){
            fillBetweenEdges(scan, AET);

            updateAET(scan, AET);

            resortAET(AET);
        }

    }

}

bool buildEdgeList(const attr_point* points, int count){

    int v1, v2;
    bool scanline_crossed = false;

    v1 = count-1;

    for(v2=0;v2<count;v2++){

        if(points[v1].coord[1] != points[v2].coord[1]){
            scanline_crossed = true;


        if(points[v1].coord[1] < points[v2].coord[1]){
            makeEdgeRec(points[v2], points[v1]);
        }
        else{
            makeEdgeRec(points[v1], points[v2]);
        }
        }
        v1 = v2;
    }

    return scanline_crossed;
}

void makeEdgeRec(attr_point upper, attr_point lower){

    float dy = upper.coord[1] - lower.coord[1];

    Edge* e = new Edge();

   if (dy != 0) {
        for(int i=0; i<ATTR_SIZE; i++){
            e->inc.coord[i] = (upper.coord[i] - lower.coord[i]) / dy;
        }
    } else {

        for(int i=0; i<ATTR_SIZE; i++){
            e->inc.coord[i] = 0;  // Set the increment to zero
        }
    }

    float factor = ceil(lower.coord[1]) - lower.coord[1];

    for(int i=0;i<ATTR_SIZE; i++){
        e->p.coord[i] = lower.coord[i] + factor*e->inc.coord[i];
    }

    e->yLast = ceil(upper.coord[1]) - 1;

    int scanline = ceil(lower.coord[1]);
    insertEdge(edge_table[scanline], e);

}
/*
void insertEdge(Edge* &list, Edge* e){

    if (list == nullptr || e->p.coord[0] < list->p.coord[0]) {
        e->next = list;
        list = e;
        return;
    }
    Edge *p, *q ;


    q=list;
    p = q->next;

    while(p!=nullptr &&(e->p.coord[0] > p->p.coord[0])){
        q=p;
        p=p->next;
    }

    e->next = q->next;
    q->next = e;

}*/

void insertEdge(Edge* &list, Edge* e) {
 /*   Edge *p = , *q = list;

     if (list == nullptr) {
        // AET is null, set list to e
        list = e;
        e->next = nullptr; // Make sure e points to nullptr
        return;
    }
    // p leads
    p = q->next;

    while (p != nullptr && (e->p.coord[0] > p->p.coord[0])) {
        // Step to the next edge
        q = p;
        p = p->next;
    }

    // Link the edge into the list after q

    e->next = q->next;
    q->next = e;
    */
    Edge *q = nullptr;
    Edge *p = list;

    while (p != nullptr && (e->p.coord[0] > p->p.coord[0])) {
        // Step to the next edge
        q = p;
        p = p->next;
    }

    if(q==nullptr){
        e->next = p;
        list = e;
    }
    else{
         e->next = q->next;
        q->next = e;
    }
}

void addActiveList(int scanline, Edge*& AET){

    Edge *p = edge_table[scanline]->next;
    Edge *q;

    while(p != nullptr){
        q = p->next;
        insertEdge(AET, p);
        p = q;
    }
    edge_table[scanline]->next = nullptr;
}

void updateAET(int scanline, Edge*& active){

    Edge *q = nullptr, *p = active;

    while(p){
        if(scanline == p->yLast){
                if(q==nullptr){
                    active = p->next;
                    delete p;
                    p =active;
                }
                else{
                    p=p->next;
                    deleteAfter(q);
                }
           // p = p->next;
            //deleteAfter(q);
        }
        else{

            for (int i = 0; i < ATTR_SIZE; ++i) {
                p->p.coord[i] += p->inc.coord[i]; // Update each coordinate separately
            }

            q=p;
            p = p->next;
        }
    }

}

void deleteAfter(Edge* q){


    Edge *p = q->next;
    q->next = p->next;
    delete p;
}

void resortAET(Edge*& active){

     Edge *sortedAET = nullptr;

     Edge *p = active;
     active = nullptr;

      while (p != nullptr) {
        Edge *q = p->next;

        insertEdge(active, p);

        p = q;
    }
   // active->next = sortedAET;
}

void fillBetweenEdges(int scanline, Edge* AET){


    Edge *p1, *p2;
    p1 = AET;

    while(p1 != nullptr && p1->next!=nullptr){
        p2 = p1->next;

        if (p2 != nullptr) {
        if(p1->p.coord[0]!=p2->p.coord[0]){
            float dx = p2->p.coord[0] - p1->p.coord[0];
            attr_point inc;
            for(int i=0;i<ATTR_SIZE;i++){
                inc.coord[i] = (p2->p.coord[i]-p1->p.coord[i])/dx;
            }
            float factor = ceil(p1->p.coord[0]) - p1->p.coord[0];
            attr_point value = p1->p;
            for (int i = 0; i < ATTR_SIZE; ++i) {
                value.coord[i] += factor * inc.coord[i];
            }

            int endx = (ceil(p2->p.coord[0]));
             while (value.coord[0] < endx) {
                // Calculate the color for the pixel and plot it
                // Here, x and z come from the current values, and y is the current scanline
               int x = static_cast<int>(value.coord[0]);
                int y = static_cast<int>(scanline);
                float z = value.coord[2];
                rd_display_pixel(x,y,z);
                // Increment the values
                for (int i = 0; i < ATTR_SIZE; ++i) {
                    value.coord[i] += inc.coord[i];
                }


               }
        }


        }
          p1=p2->next;
    }
}
int poly_clip(int n_vertex, attr_point vertex_list[], attr_point clipped_list[]) {
/*
     for (int i = 0; i < n_vertex; ++i) {
        clipped_list[i] = vertex_list[i];
    }
    return n_vertex;
*/
    int count=0;
     attr_point first[6];
     bool flag[6] = {false};
     attr_point last[6];

    for(int i=0;i<n_vertex;i++){
        clip_a_point(vertex_list[i], Boundary::LEFT, first, last, flag, clipped_list, count);
    }

    clip_last_point(first, last, flag, clipped_list, count);

    return count;
}

void clip_a_point(const attr_point& p, Boundary b, attr_point (&first)[], attr_point (&last)[], bool (&flag)[], attr_point clipped_list[], int& count){

    if(!flag[static_cast<int>(b)]){
        first[static_cast<int>(b)] = p;
        flag[static_cast<int>(b)] = true;
    }else{
        if(cross(p, last[static_cast<int>(b)], b)){

            attr_point ipt = intersect(p, last[static_cast<int>(b)], b);

            if(b != Boundary::BACK){
                clip_a_point(ipt,  static_cast<Boundary>(static_cast<int>(b) + 1), first,last,flag,clipped_list,count);
            }
            else{
                clipped_list[count++] = ipt;
            }
        }
    }

    last[static_cast<int>(b)] = p;

    if(inside(p,b)){
        if(b != Boundary::BACK){
            clip_a_point(p, static_cast<Boundary>(static_cast<int>(b)+1), first, last, flag, clipped_list, count);
        }
        else{
            clipped_list[count++] = p;
        }
    }
}

void clip_last_point(attr_point (&first)[], attr_point (&last)[], bool (&flag)[], attr_point clipped_list[], int& count){

    for(int b=0;b<6;b++){
        if(flag[b]){
            if(cross(first[b],last[b],static_cast<Boundary>(b))){
                attr_point ipt = intersect(last[b], first[b], static_cast<Boundary>(b));

                if(b != Boundary::BACK){
                clip_a_point(ipt, static_cast<Boundary>(static_cast<int>(b)+1), first, last, flag, clipped_list, count);
                }
                else{
                    clipped_list[count++] = ipt;
                }
            }
        }
    }

}

bool cross(const attr_point& p1, const attr_point& p2, Boundary b){

    bool inside_p1 = inside(p1, b);
    bool inside_p2 = inside(p2, b);

    return inside_p1 != inside_p2;
}

bool inside(const attr_point& p, Boundary b){
    float bc;

    switch(b){
    case Boundary::LEFT:
        bc = p.coord[0];
        break;
    case Boundary::RIGHT:
        bc = p.coord[3]-p.coord[0];
        break;
    case Boundary::BOTTOM:
        bc = p.coord[1];
        break;
    case Boundary::TOP:
        bc = p.coord[3]-p.coord[1];
        break;
    case Boundary::FRONT:
        bc = p.coord[2];
        break;
    case Boundary::BACK:
        bc = p.coord[3]-p.coord[2];
        break;
    default:
        return false;
    }

    return (bc >= 0.0f );

}

attr_point intersect(const attr_point& p1, const attr_point& p2, Boundary b){

    attr_point point = {0.0f,0.0f,0.0f,1.0f};
    float alpha;

    switch(b){
    case Boundary::LEFT:
        alpha = p1.coord[0] / (p1.coord[0]-p2.coord[0]);
        break;
    case Boundary::RIGHT:
        alpha = (p1.coord[3]-p1.coord[0]) / ((p1.coord[3]-p1.coord[0])- (p2.coord[3]-p2.coord[0]));
        break;
    case Boundary::BOTTOM:
        alpha = p1.coord[1] / (p1.coord[1]-p2.coord[1]);
        break;
    case Boundary::TOP:
        alpha = (p1.coord[3]-p1.coord[1]) / ((p1.coord[3]-p1.coord[1])- (p2.coord[3]-p2.coord[1]));
        break;
    case Boundary::FRONT:
        alpha = p1.coord[2] / (p1.coord[2]-p2.coord[2]);
        break;
    case Boundary::BACK:
        alpha = (p1.coord[3]-p1.coord[2]) / ((p1.coord[3]-p1.coord[2])- (p2.coord[3]-p2.coord[2]));
        break;
    default:
        return point;
    }

    for(int i=0;i<ATTR_SIZE;i++){
        point.coord[i] = p1.coord[i] + alpha*(p2.coord[i]-p1.coord[i]);
    }

    return point;
}
