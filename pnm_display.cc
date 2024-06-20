#include <iostream>
#include <fstream>
#include "pnm_display.h"
#include "rd_error.h"
#include "rd_display.h"
#include <string>

#define RGB_VALUE 3

using std::string;
using namespace std;

float **image= nullptr;
int pnm_frame_number = 0;
float background_color[3] = {0.0f,0.0f,0.0f};

int pnm_init_display(){

    if(display_xSize <= 0 || display_ySize <=0){

        return RD_INPUT_UNINITIALIZED_DISPLAY;
    }

    image = new float*[display_ySize];

    for(int i=0;i<display_ySize;i++){

        image[i] = new float[display_xSize*RGB_VALUE];

    }

    return RD_OK;
}

int pnm_end_display(){

    for (int i=0; i< display_ySize;i++){
        delete[] image[i];
    }
    delete[] image;

    return RD_OK;
}

int pnm_init_frame(int frame_no){

    pnm_frame_number = frame_no;

  /*  for(int i=0;i<display_xSize;i++){
        for(int j=0;j<display_ySize;j++){
            image[j][i] = 0.0f;
        }
    }
*/
    pnm_clear();

    return RD_OK;
}

int pnm_end_frame(){

    string file_name = display_name+std::to_string(pnm_frame_number)+".ppm";

    ofstream outputFile(file_name);
    if(!outputFile){
        return RD_INPUT_DISPLAY_INITIALIZATION_ERROR;
    }

    //header
    outputFile << "P3" << endl;
    outputFile << display_xSize << " " << display_ySize << endl;
    outputFile << 255 << endl;

    for(int i =0; i< display_ySize; i++){

        for(int j=0;j<display_xSize;j++){

            outputFile << static_cast<int>(image[i][j * RGB_VALUE] * 255) << " ";
            outputFile << static_cast<int>(image[i][j * RGB_VALUE +1] * 255) << " ";
            outputFile << static_cast<int>(image[i][j * RGB_VALUE +2] * 255) << " ";

        }
        outputFile << endl;
    }

    outputFile.close();
    return RD_OK;
}

int pnm_write_pixel(int x, int y, const float rgb []){

    if(x <0 || x>= display_xSize || y <0 || y>= display_ySize){
        return RD_OK;
    }

    image[y][x * RGB_VALUE] = rgb[0] > 1.0f ? 1.0f : (rgb[0] < 0.0f ? 0.0f : rgb[0]);
    image[y][x * RGB_VALUE + 1] = rgb[1] > 1.0f ? 1.0f : (rgb[1] < 0.0f ? 0.0f : rgb[1]);
    image[y][x * RGB_VALUE + 2] = rgb[2] > 1.0f ? 1.0f : (rgb[2] < 0.0f ? 0.0f : rgb[2]);

    return RD_OK;
}

int pnm_read_pixel(int x, int y, float rgb []){

    if(x <0 || x>= display_xSize || y <0 || y>= display_ySize){
        return RD_OK;
    }

    rgb[0] = image[y][x * RGB_VALUE] > 1.0f ? 1.0f : (image[y][x * RGB_VALUE] < 0.0f ? 0.0f : image[y][x * RGB_VALUE]);
    rgb[1] = image[y][x * RGB_VALUE + 1] > 1.0f ? 1.0f : (image[y][x * RGB_VALUE + 1] < 0.0f ? 0.0f : image[y][x * RGB_VALUE + 1]);
    rgb[2] = image[y][x * RGB_VALUE + 2] > 1.0f ? 1.0f : (image[y][x * RGB_VALUE + 2] < 0.0f ? 0.0f : image[y][x * RGB_VALUE + 2]);


    return RD_OK;
}

int pnm_set_background(const float rgb []){

    background_color[0] = rgb[0] > 1.0f ? 1.0f : (rgb[0] < 0.0f ? 0.0f : rgb[0]);
    background_color[1] = rgb[1] > 1.0f ? 1.0f : (rgb[1] < 0.0f ? 0.0f : rgb[1]);
    background_color[2] = rgb[2] > 1.0f ? 1.0f : (rgb[2] < 0.0f ? 0.0f : rgb[2]);

    return RD_OK;
}

int pnm_clear(){

    for (int i = 0; i < display_ySize; i++) {
        for (int j = 0; j < display_xSize * RGB_VALUE; j++) {

            image[i][j] = background_color[j % RGB_VALUE];
        }
    }
    return RD_OK;
}
