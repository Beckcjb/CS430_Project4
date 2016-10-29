// Charles Beck
// PPM Writer
#include <stdio.h>
#include <stdlib.h>

typedef struct PixelColor {
    unsigned char r, g, b;
} PixelColor;


typedef struct Pixmap
{
    int width, height, magicNumber, color;
    PixelColor *image;
}Pixmap;

int ppmWriter(Pixmap *buffer, char *outputFileName, int size, int type)
{
    FILE *printFile;
    int i, j, numPix;
    char comment[] = {"#Creator: Charles Beck"};


    printFile = fopen(outputFileName, "w");
    if (!printFile)
    {
        fprintf(stderr,"Erroe: Can't open the file for writing");
        fclose(printFile);
        exit(1);
    }
    else
    {
        fprintf(printFile, "P%d\n%s\n%d %d\n%d\n", type, comment, buffer->width, buffer->height, 255);
     
											// Print out to the outfile in P3 format
        if(type == 3)
        {
            for(i = 0; i < (buffer->height); i++)
            {
                for(j = 0; j < (buffer->width); j++)
                {
                    fprintf(printFile, "%d ", buffer->image[i * buffer->width *3+3*j].r);
                    fprintf(printFile, "%d ", buffer->image[i * buffer->width *3+3*j+1].g);
                    fprintf(printFile, "%d\n", buffer->image[i * buffer->width *3+ 3*j+2].b);
                }
            }
        }
    }
    fclose(printFile);
    return 0;
}
