import numpy as np
import matplotlib.pyplot as plt
import os

### Isn't this square matrix and other one sparse matrix???

def readSquareMatrix(filename, total_length, resolution):
    infile = open(filename).readlines()
    HiC = np.zeros((total_length,total_length)).astype(np.int16)
    percentage_finish = 0
    for i in range(0, len(infile)):
        if (i %  (len(infile) / 10)== 0):
            print ('finish ', percentage_finish, '%')
            percentage_finish += 10
        nums = infile[i].split()
        x = int(int(nums[0])/resolution)
        y = int(int(nums[1])/resolution)
        val = int(float(nums[2]))

        HiC[x][y] = val
        HiC[y][x] = val
    return HiC


### this function return result and index. result is array with shape (N,40,40) including N 40*40 frames from chrN. index is array with shape (N,4) including
### corresponding information for each frame (tag, chr_number, i, j). i and j are indices of cell[0,0] in the frame.
def divide(HiCfile, chrN, input_resolution):
    subImage_size = 40
    step = 25
    chrs_length = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566]
    result = []
    index = []
    matrix_name = HiCfile + '_npy_form_tmp.npy'
    if os.path.exists(matrix_name):
        print ('loading ', matrix_name)
        HiCsample = np.load(matrix_name)
    else:
        #print (matrix_name, 'not exist, creating ', HiCfile)
        print(chrN)
        HiCsample = readSquareMatrix(HiCfile, (int(chrs_length[chrN-1]/input_resolution) + 1), input_resolution)
        np.save(matrix_name, HiCsample)

    ### need to think more about that
    total_loci = int(chrs_length[chrN-1] / input_resolution)
    for i in range(0, total_loci - subImage_size, step):
        for j in range(np.maximum(0,i-200), np.minimum(i+201,total_loci-subImage_size), step): # not better to have step here too???
            subImage = HiCsample[i:i+subImage_size, j:j+subImage_size]
            # print("i: ", i, ", j: ", j)
            result.append(subImage)
            tag = 'test'
            index.append((tag, chrN, i, j))
            ### testing result2 is made correctly
            #if i == 925 and j == 941:
            #    rr = len(result)
    print(len(result))
    #result = np.array(result)
    #index = np.array(index)
    ### testing result is made correctly, it should be = 786 for index i = 925 and j = 941
    # print(result[rr-1,16,0])
    # print(index[rr-1,])
    return result, index


def divide2(HiC_mat):
    subImage_size = 40
    step = 25
    total_loci = HiC_mat.shape[0]
    result = []
    index = []
    for i in range(0, total_loci - subImage_size, step):
        for j in range(np.maximum(0,i-200), np.minimum(i+200, total_loci - subImage_size), step):
            subImage = HiC_mat[i:i + subImage_size, j:j + subImage_size]
            result.append(subImage)
            index.append((i,j))
    return result, index


if __name__ == "__main__":
    main()
