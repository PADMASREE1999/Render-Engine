#ifndef RD_ERROR_H
#define RD_ERROR_H

#include <string>
using std::string;

/************** Some error codes ************************************/

#define RD_OK                                     0
#define RD_INPUT_EXPECTED_EOF                    -1
/* File ended, but in a good place */

#define RD_INPUT_UNINITIALIZED_DISPLAY           -2
#define RD_INPUT_UNAVAILABLE_DISPLAY             -3
#define RD_INPUT_INCOMPATIBLE_DISPLAY            -4
#define RD_INPUT_DISPLAY_INITIALIZATION_ERROR    -5
/* Can't get a full color display */

#define RD_INPUT_UNEXPECTED_EOF                  -6
#define RD_INPUT_FILE_NOT_FOUND                  -7
#define RD_INPUT_UNKNOWN_TOKEN                   -8
/*  Completely unrecognized token */

#define RD_INPUT_ILLEGAL_TOKEN                   -9
/*  Recognized token, but in the wrong place */
#define RD_INPUT_ILLEGAL_NAN                    -10
/*  NaN in unacceptable location */

#define RD_INPUT_RUNAWAY_STRING                 -11
#define RD_INPUT_UNKNOWN_AXIS                   -12
#define RD_INPUT_OBJ_LIST_OVERFLOW              -13
#define RD_INPUT_CIRCULAR_INSTANCE              -14
#define RD_INPUT_UNKNOWN_VERTEX_TYPE            -15
#define RD_INPUT_ILLEGAL_VERTEX_INDEX           -16
#define RD_INPUT_UNKNOWN_OBJECT_INSTANCE        -17
#define RD_INPUT_DUPLICATE_OBJ_LABEL            -18
#define RD_INPUT_ILLEGAL_FLAG_VALUE             -19
#define RD_INPUT_INSUFFICIENT_DISPLAY_MEMORY    -20
#define RD_INPUT_INSUFFICIENT_ZBUFFER_MEMORY    -21
#define RD_INPUT_INSUFFICIENT_TEXTURE_MEMORY    -22
#define RD_INPUT_INTERLEAVED_BLOCK              -23
/* The end of this block is not the type expected */

#define RD_INPUT_UNKNOWN_DISPLAY_TYPE           -24
#define RD_INPUT_UNKNOWN_DISPLAY_MODE           -25

#define RD_INPUT_UNKNOWN_SURFACE_TYPE           -26
#define RD_INPUT_LIGHT_LIST_OVERFLOW            -27

#define RD_INPUT_UNKNOWN_CURVE_TYPE             -28
#define RD_INPUT_UNKNOWN_PATCH_TYPE             -29

#define RD_INPUT_TRANSFORM_STACK_OVERFLOW       -30
#define RD_INPUT_TRANSFORM_STACK_UNDERFLOW      -31

#define RD_INPUT_MAP_LIST_OVERFLOW              -32
#define RD_INPUT_MAP_FILE_NOT_FOUND             -33
#define RD_INPUT_MAP_BAD_FORMAT                 -34
#define RD_INPUT_MAP_ILLEGAL_LABEL              -35
#define RD_INPUT_MAP_NOT_FOUND                  -36
#define RD_INPUT_MAP_UNKNOWN_MAP_TYPE           -37
#define RD_INPUT_MAP_BAD_INTERP                 -38
#define RD_INPUT_MAP_BAD_BORDER                 -39

#define RD_INPUT_UNKNOWN_ENGINE                 -40

#define RD_INPUT_NON_INTEGER_VALUE              -41
#define RD_INPUT_EXPECTED_VERTEX_INDEX          -42
#define RD_INPUT_EXPECTED_VERTEX_COORD          -43
#define RD_INPUT_EXPECTED_REAL                  -44

#define RD_INPUT_UNKNOWN_OPTION                 -45

#define RD_INPUT_BAD_BSPLINE_DEGREE             -46

#define RD_INPUT_BAD_OBJECT_BLOCK_PARAMETER     -47
#define RD_INPUT_BAD_OBJECT_INSTANCE_PARAMETER  -48

#if defined(__cplusplus)
extern "C"
{
#endif

  void rd_print_error(int err, const string & filename);

#if defined (__cplusplus)
}
#endif

#endif /* RD_ERROR_H */
