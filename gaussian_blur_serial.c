
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846


int main(int argc, char * argv[])
{   
    /**To ensure there are three arguments: <input_pgm>, <output_pgm> and <sigma>*/
    if (argc != 4) {
         fprintf(stderr, "Usage: ./gaussian_blur_serial <input_pgm> <output_pgm> <sigma>\n");
    }

    /**Parse the command line to recieve the input file, output file name and the sigma value*/
    const char* input_file = argv[1]; //input binary PGM file
    const char* output_file = argv[2]; //output binary PGM file
    float sigma = atof(argv[3]); //sigma value

    size_t width, height;
    unsigned char max_value;

    /**Open the input binary PGM file*/
    FILE* file = fopen(input_file, "rb");
    if (file == NULL) {
        perror("Cannot open file");
        exit(1);
    }

    //read through the input PGM file
    if (fscanf(file, "P5\n%ld %ld\n%hhd\n", &width, &height, &max_value) != 3) {
        fprintf(stderr, "Error reading image header\n");
        fclose(file);
        return EXIT_FAILURE;
    }

    unsigned char * image = malloc(width*height); // image binary data
    unsigned char * blurred_image = malloc(width*height); // blurred image
    /** Push image into buffer */
    int seek_dist = 3 + (int)log10(height) + (int)log10(width) + 4 + 4; // |firstLine| + |secondLine| + |thirdLine|
    fseek(file, seek_dist, SEEK_SET);
    if (fread(image, 1, height * width, file) != height * width) {
        fprintf(stderr, "Error reading image data\n");
        free(image);
        fclose(file);
        return 1;
    }


    /**Validate the bounds of sigma*/
    if (sigma <= 0) {
        fprintf(stderr, "Sigma value must be greater than 0\n");
        exit(1);
    }
    
    /**Create the gaussian kernel matrix*/
    uint32_t kernel_order = ceil(sigma * 6); //order of the kernel matrix
    uint32_t half_order = kernel_order / 2;

     /**Validate the bounds of the kernel matrix*/
    if (kernel_order > width || kernel_order > height) {
        fprintf(stderr, "The kernel matrix should not be bigger than the input image size \n");
        return EXIT_FAILURE;
    }

    if (kernel_order % 2 == 0) kernel_order += 1; //make it odd to account for edges

    //float** kernel_matrix = create_kernel_matrix(kernel_order, sigma);
    float kernel_matrix[kernel_order][kernel_order];
    float normalizer = 1 / (2 * PI * sigma * sigma);
    float sum = 0.0;

    //fills in the gaussian kernel matrix
    for (size_t x = 0; x < kernel_order; x++) {
        for (size_t y = 0; y < kernel_order; y++) {
            uint32_t x_offset = x - half_order;
            uint32_t y_offset = y - half_order;
            kernel_matrix[x][y] = normalizer * exp(-1.0 * (x_offset * x_offset + y_offset * y_offset) / (2.0 * sigma * sigma));
            sum += kernel_matrix[x][y];
        }
    }

    /** normalizes the gaussian kernel matrix to the total sum */
    for (size_t i = 0; i < kernel_order; i++) {
        for (size_t j = 0; j < kernel_order; j++) {
            kernel_matrix[i][j] /= sum;
        }
    }

    /** apply convolution to all pixels in image */
      size_t dc = (kernel_order - 1) / 2;
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            float conv = 0.0;
            float weight_sum = 0.0;
            for (size_t ki = 0; ki < kernel_order; ki++) {
                for (size_t kj = 0; kj < kernel_order; kj++) {
                    int ii = i + ki - dc;
                    int jj = j + kj - dc;

                    // Handle boundaries by clamping to the nearest edge pixel
                    if (ii < 0) ii = 0;
                    if (ii >= height) ii = height - 1;
                    if (jj < 0) jj = 0;
                    if (jj >= width) jj = width - 1;

                    conv += image[ii * width + jj] * kernel_matrix[ki][kj];
                    weight_sum += kernel_matrix[ki][kj];
                }
            }
            // Normalize the convolution result to prevent overflow
            if (weight_sum > 0) {
                conv /= weight_sum;
            }
            // Ensure the result is within the valid range for unsigned char
            blurred_image[i * width + j] = (unsigned char)(conv);
        }
    }

    /** output file */
    FILE * out = fopen(output_file, "wb");
    if(!out){
        fprintf(stderr, "Error with opening file!\n");
        exit(1);
    }
    fprintf(out, "P5\n%ld %ld\n%d\n", width, height, max_value); 

    fwrite(blurred_image, sizeof(char), width*height, out);

    fclose(out);
    fclose(file);
    free(image);
    free(blurred_image);



    return 0;

} 