#ifndef COLOR_H
#define COLOR_H

#include <vector>
#include <iostream>

/*
   Return a RGB colour value given a scalar v in the range [vmin,vmax]
   In this case each colour component ranges from 0 (no contribution) to
   1 (fully saturated), modifications for other ranges is trivial.
   The colour is clipped at the end of the scales if v is outside
   the range [vmin,vmax]
*/

typedef struct {
    double a,b,c;
    inline void print() { std::cout << a << " " << b << " " << c << std::endl; }
} color;

color get_rgb_color(double v,double vmin,double vmax);
void get_rgb_palette(std::vector<color> &palette, int palette_size);
void get_hsv_palette(std::vector<color> &palette, int palette_size);
void print_palette(std::vector<color> &palette);
color rgb_2_hsv(color &rgb);


#endif // COLOR_H
