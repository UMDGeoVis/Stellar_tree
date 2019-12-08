#include "color.h"

color get_rgb_color(double v,double vmin,double vmax)
{
   color c = {1.0,1.0,1.0}; // white
   double dv;

   if (v < vmin)
      v = vmin;
   if (v > vmax)
      v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) {
      c.a = 0;
      c.b = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5 * dv)) {
      c.a = 0;
      c.c = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   } else if (v < (vmin + 0.75 * dv)) {
      c.a = 4 * (v - vmin - 0.5 * dv) / dv;
      c.c = 0;
   } else {
      c.b = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
      c.c = 0;
   }

   return c;
}

void get_rgb_palette(std::vector<color> &palette, int palette_size)
{
    double palette_min = 0, palette_max = 1.0;
    double incr_palette = 1.0/palette_size;

//    std::cout << "incr_palette: " << incr_palette << std::endl;

    for(int i=0; i<palette_size; i++)
    {
        double value = palette_min + i * incr_palette;
        color c = get_rgb_color(value,palette_min,palette_max);
        palette.push_back(c);
    }
}

void get_hsv_palette(std::vector<color> &palette, int palette_size)
{
    double palette_min = 0, palette_max = 1.0;
    double incr_palette = 1.0/palette_size;

//    std::cout << "incr_palette: " << incr_palette << std::endl;

    for(int i=0; i<palette_size; i++)
    {
        double value = palette_min + i * incr_palette;
        color c = get_rgb_color(value,palette_min,palette_max);
//        std::cout<<"rgb "; c.print();
        color h = rgb_2_hsv(c);
//        std::cout<<"hsv "; h.print();
        palette.push_back(/*rgb_2_hsv(c)*/h);
    }
}

void print_palette(std::vector<color> &palette)
{
    for(unsigned i=0; i<palette.size(); i++)
    {
        color &c = palette[i];
        std::cout<<c.a<<" "<<c.b<<" "<<c.c<<std::endl;
    }
    int a; std::cin >> a;
}

color rgb_2_hsv(color &rgb)
{
    color hsv;
    double rgb_max = std::max(rgb.a, std::max(rgb.b, rgb.c));
    double rgb_min = std::min(rgb.a, std::min(rgb.b, rgb.c));
    double delta = rgb_max - rgb_min;
    hsv.b = delta / (rgb_max + 1e-20d);
    hsv.c = rgb_max;

    double hue;
    if (rgb.a == rgb_max)
        hue = (rgb.b - rgb.c) / (delta + 1e-20d);
    else if (rgb.b == rgb_max)
        hue = 2 + (rgb.c - rgb.a) / (delta + 1e-20d);
    else
        hue = 4 + (rgb.a - rgb.b) / (delta + 1e-20d);
    if (hue < 0)
        hue += 6.d;
    hsv.a = hue * (1.d / 6.d);
    return hsv;
}
