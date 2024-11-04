#ifndef IMAGELOADER_H
#define IMAGELOADER_H

#include <string>
#include <vector>

// Source courtesy of J. Manson
// http://josiahmanson.com/prose/optimize_ppm/

namespace ppmLoader {
    using namespace std;

    /**
     * Skips comments in the PPM file.
     * @param f The input file stream.
     */
    void eat_comment(ifstream &f);

    /**
     * Structure representing an RGB color.
     */
    struct RGB {
        unsigned char r;    ///< Red component.
        unsigned char g;    ///< Green component.
        unsigned char b;    ///< Blue component.
    };

    /**
     * Structure representing an RGB image.
     */
    struct ImageRGB {
        int w;      ///< Width of the image.
        int h;      ///< Height of the image.
        vector<RGB> data; ///< Pixel data of the image.
    };

    /**
     * Loads a PPM image into an ImageRGB structure.
     * @param img The ImageRGB structure to load the image into.
     * @param name The name of the PPM file.
     */
    void load_ppm(ImageRGB &img, const string &name);

    /**
     * Enum representing the format of the loaded image.
     */
    enum loadedFormat {
        rgb,    ///< RGB format.
        rbg     ///< RBG format.
    };

    /**
     * Loads a PPM image into a pixel array.
     * @param pixels The pixel array to load the image into.
     * @param w The width of the image.
     * @param h The height of the image.
     * @param name The name of the PPM file.
     * @param format The format of the loaded image (default is RGB).
     */
    void load_ppm(
        unsigned char *&pixels, unsigned int &w, unsigned int &h,
        const string &name, loadedFormat format = rgb
    );
} // namespace ppmLoader

#endif // IMAGELOADER_H