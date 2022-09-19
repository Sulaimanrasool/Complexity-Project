import numpy as np
import matplotlib.pyplot as plt
import random as rand

"This File is the code for the Oslo Model"


def oslo_algorithm(L,g):
    """

    :param L: The number of sites
    :param g : grains added
    :return steady_heights_array: 2D array of heights , L * g , after steady state has been reached
    :return avalanche_array_steady : Array of avalanche sizes after steady state has been reached
    """
    thresh_slope = np.random.randint(low = 1,high =  3, size  = L)  # Random threshold values
    heights = np.zeros((1, L + 1))        # Create an array of heights (zeros) for all sites up till l + 1(always zero)
    slope = -np.diff(heights, n = 1, axis= -1)  # Slope of each site
    avalanche_array_steady = []
    avalanche_array_tran = []
    steady_heights_array = np.empty((0,L + 1), int)
    tran_heights_array = np.empty((0, L + 1), int)
    slope_array = np.empty((0, L), int)

    grains_left = g
    steady_state = 0
    while grains_left > 0 :

        heights[0][0] = heights[0][0] + 1             # add a grain to the first site
        slope[0][0] = slope[0][0] + 1                 # add 1 to the slope
        relaxation = 0
        grains_left -= 1
        while (slope > thresh_slope).any():           # while there is a slope with a value greater than threshold

            for i in range(L):

                if slope[0][i] > thresh_slope[i]:
                    if i == 0 :
                        slope[0][0] -= 2    # decrease the slope of first site
                        slope[0][1] += 1    # increase the slope of the second site
                        heights[0][0] -= 1  # decrease the height of site 1
                        heights[0][1] += 1  # increase the height of site 2
                        thresh_slope[0] = rand.randint(1,2)
                        relaxation += 1     # add one to the relaxation


                    elif 0 < i <= L-2 :      # for all sites not the first or last
                        slope[0][i] -= 2    # decrease slope of site i
                        slope[0][i+1] += 1      # increase slope of site i + 1
                        slope[0][i-1] += 1      # increase the slope of site i - 1
                        heights[0][i] -= 1      # decrease the height of site i
                        heights[0][i+1] += 1    # increase the height of site i + 1
                        thresh_slope[i] = rand.randint(1, 2)
                        relaxation += 1         # add one to the relaxation

                    if i == L-1 :               # for the last site
                        slope[0][i] -= 1        # minus 1 as the last pieces falls onto a zero area
                        slope[0][i-1] += 1     # site before last site
                        heights[0][i] -= 1      # take away 1 from the height
                        thresh_slope[i] = rand.randint(1, 2)
                        relaxation += 1
                        if steady_state == 0:   # tell me when steady state has been reached but not after
                            print("Steady state has been reached, grains required to reach steady state is : ",  g-grains_left)
                            avalanche_array_tran.append(relaxation)
                            tran_heights_array = np.append(tran_heights_array, np.array(heights), axis=0)
                            steady_state += 1
        if steady_state == 0:
            tran_heights_array = np.append(tran_heights_array, np.array(heights), axis = 0)
            avalanche_array_tran.append(relaxation)
        else:
            if grains_left > 0:
                steady_heights_array = np.append(steady_heights_array, np.array(heights), axis = 0)
                avalanche_array_steady.append(relaxation)
                slope_array = np.append(slope_array, slope, axis = 0)


    return steady_heights_array, avalanche_array_steady, tran_heights_array, slope_array

def oslo_algorithm_P(L,g, P):
    """

    :param L: The number of sites
    :param g : grains added
    :return steady_heights_array: 2D array of heights , L * g , after steady state has been reached
    :return avalanche_array_steady : Array of avalanche sizes after steady state has been reached
    """
    # np.random.seed(seed = 1)   # Seed control
    thresh_slope = np.full((1,L), P, dtype=int) # Random threshold values
    # plt.hist(thresh_slope)
    # plt.show()
    heights = np.zeros((1, L + 1))        # Create an array of heights (zeros) for all sites up till l + 1(always zero)
    slope = -np.diff(heights, n = 1, axis= -1)  # Slope of each site
    avalanche_array_steady = []
    avalanche_array_tran = []
    steady_heights_array = np.empty((0,L + 1), int)
    tran_heights_array = np.empty((0, L + 1), int)

    grains_left = g
    steady_state = 0
    while grains_left > 0 :

        heights[0][0] = heights[0][0] + 1             # add a grain to the first site
        slope[0][0] = slope[0][0] + 1                 # add 1 to the slope
        relaxation = 0
        grains_left -= 1
        while (slope > thresh_slope).any():           # while there is a slope with a value greater than threshold

            for i in range(L):

                if slope[0][i] > thresh_slope[0][i]:
                    if i == 0 :
                        slope[0][0] -= 2    # decrease the slope of first site
                        slope[0][1] += 1    # increase the slope of the second site
                        heights[0][0] -= 1  # decrease the height of site 1
                        heights[0][1] += 1  # increase the height of site 2
                        thresh_slope[0][1] = P # new threshold value for slope
                        relaxation += 1     # add one to the relaxation

                    elif 0 < i <= L-2 :      # for all sites not the first or last
                        slope[0][i] -= 2    # decrease slope of site i
                        slope[0][i+1] += 1      # increase slope of site i + 1
                        slope[0][i-1] += 1      # increase the slope of site i - 1
                        heights[0][i] -= 1      # decrease the height of site i
                        heights[0][i+1] += 1    # increase the height of site i + 1
                        thresh_slope[0][i] = P  # new threshold value for slope
                        relaxation += 1         # add one to the relaxation

                    if i == L-1 :               # for the last site

                        slope[0][i] -= 1        # minus 1 as the last pieces falls onto a zero area
                        slope[0][i -1] += 1     # site before last site
                        heights[0][i] -= 1      # take away 1 from the height
                        thresh_slope[0][i] = P  # new threshold value for slope
                        relaxation += 1
                        if steady_state == 0:   # tell me when steady state has been reached but not after
                            print("Steady state has been reached, grains required to reach steady state is : ",  g-grains_left)
                            avalanche_array_tran.append(relaxation)
                            tran_heights_array = np.append(tran_heights_array, np.array(heights), axis=0)
                            steady_state += 1



        if steady_state == 0:
            tran_heights_array = np.append(tran_heights_array, np.array(heights), axis = 0)
            avalanche_array_tran.append(relaxation)
        else:
            if grains_left > 0:
                steady_heights_array = np.append(steady_heights_array, np.array(heights), axis = 0)
                avalanche_array_steady.append(relaxation)

    return steady_heights_array, avalanche_array_steady, tran_heights_array

def oslo_algorithm_ns(L,g):
    """
    No slope algorithm

    :param L: The number of sites
    :param g : grains added
    :return steady_heights_array: 2D array of heights , L * g , after steady state has been reached
    :return avalanche_array_steady : Array of avalanche sizes after steady state has been reached
    """
    # np.random.seed(seed = 1)   # Seed control
    thresh_slope = np.random.randint(low = 1,high =  3, size  = L)  # Random threshold values
    # plt.hist(thresh_slope)
    # plt.show()
    heights = np.zeros((1, L + 1))        # Create an array of heights (zeros) for all sites up till l + 1(always zero)
    slope = -np.diff(heights, n = 1, axis= -1)  # Slope of each site
    avalanche_array_steady = []
    avalanche_array_tran = []
    total_array = np.zeros((1,g))


    grains_left = g
    steady_state = 0
    val = 0
    while grains_left > 0 :

        heights[0][0] = heights[0][0] + 1             # add a grain to the first site
        slope[0][0] = slope[0][0] + 1                 # add 1 to the slope
        relaxation = 0
        grains_left -= 1
        while (slope > thresh_slope).any():           # while there is a slope with a value greater than threshold

            for i in range(L):

                if slope[0][i] > thresh_slope[i]:
                    if i == 0 :
                        slope[0][0] -= 2    # decrease the slope of first site
                        slope[0][1] += 1    # increase the slope of the second site
                        heights[0][0] -= 1  # decrease the height of site 1
                        heights[0][1] += 1  # increase the height of site 2
                        thresh_slope[0] = rand.randint(1,2)
                        relaxation += 1     # add one to the relaxation

                    elif 0 < i <= L-2 :      # for all sites not the first or last
                        slope[0][i] -= 2    # decrease slope of site i
                        slope[0][i+1] += 1      # increase slope of site i + 1
                        slope[0][i-1] += 1      # increase the slope of site i - 1
                        heights[0][i] -= 1      # decrease the height of site i
                        heights[0][i+1] += 1    # increase the height of site i + 1
                        thresh_slope[i] = rand.randint(1, 2)
                        relaxation += 1         # add one to the relaxation

                    if i == L-1 :               # for the last site
                        slope[0][i] -= 1        # minus 1 as the last pieces falls onto a zero area
                        slope[0][i-1] += 1     # site before last site
                        heights[0][i] -= 1      # take away 1 from the height
                        thresh_slope[i] = rand.randint(1, 2)
                        relaxation += 1
                        if steady_state == 0:   # tell me when steady state has been reached but not after
                            print("Steady state has been reached, grains required to reach steady state is : ",  g-grains_left)
                            avalanche_array_tran.append(relaxation)
                            total_array[0][val] = heights[0][0]
                            tran_heights_array = total_array[0][0:val]

                            steady_state += 1

                            sted = val
                            val += 1

        if steady_state == 0:
            total_array[0][val] = heights[0][0]
            val += 1
            avalanche_array_tran.append(relaxation)
        else:
            if grains_left > 0:
                total_array[0][val] = heights[0][0]
                avalanche_array_steady.append(relaxation)
                val += 1

    steady_heights_array = total_array[0][sted:]
    return steady_heights_array, avalanche_array_steady, tran_heights_array,  total_array

def oslo_algorithm_cross(L,g):
    """
    cross over time algorithm

    :param L: The number of sites
    :param g : grains added
    :return steady_heights_array: 2D array of heights , L * g , after steady state has been reached
    :return avalanche_array_steady : Array of avalanche sizes after steady state has been reached
    """
    # np.random.seed(seed = 1)   # Seed control
    thresh_slope = np.random.randint(low = 1,high =  3, size  = L)  # Random threshold values
    # plt.hist(thresh_slope)
    # plt.show()
    heights = np.zeros((1, L + 1))        # Create an array of heights (zeros) for all sites up till l + 1(always zero)
    slope = -np.diff(heights, n = 1, axis= -1)  # Slope of each site
    avalanche_array_steady = []
    avalanche_array_tran = []
    total_array = np.zeros((1,g))


    grains_left = g
    steady_state = 0
    val = 0
    while grains_left > 0 :

        heights[0][0] = heights[0][0] + 1             # add a grain to the first site
        slope[0][0] = slope[0][0] + 1                 # add 1 to the slope
        relaxation = 0
        grains_left -= 1
        while (slope > thresh_slope).any():           # while there is a slope with a value greater than threshold

            for i in range(L):

                if slope[0][i] > thresh_slope[i]:
                    if i == 0 :
                        slope[0][0] -= 2    # decrease the slope of first site
                        slope[0][1] += 1    # increase the slope of the second site
                        heights[0][0] -= 1  # decrease the height of site 1
                        heights[0][1] += 1  # increase the height of site 2
                        thresh_slope[0] = rand.randint(1,2)
                        relaxation += 1     # add one to the relaxation

                    elif 0 < i <= L-2 :      # for all sites not the first or last
                        slope[0][i] -= 2    # decrease slope of site i
                        slope[0][i+1] += 1      # increase slope of site i + 1
                        slope[0][i-1] += 1      # increase the slope of site i - 1
                        heights[0][i] -= 1      # decrease the height of site i
                        heights[0][i+1] += 1    # increase the height of site i + 1
                        thresh_slope[i] = rand.randint(1, 2)
                        relaxation += 1         # add one to the relaxation

                    if i == L-1 :               # for the last site
                        slope[0][i] -= 1        # minus 1 as the last pieces falls onto a zero area
                        slope[0][i-1] += 1     # site before last site
                        heights[0][i] -= 1      # take away 1 from the height
                        thresh_slope[i] = rand.randint(1, 2)
                        relaxation += 1
                        if steady_state == 0:   # tell me when steady state has been reached but not after
                            print("Steady state has been reached, grains required to reach steady state is : ",  g-grains_left)

                            avalanche_array_tran.append(relaxation)
                            total_array[0][val] = heights[0][0]
                            tran_heights_array = total_array[0][0:val]

                            steady_state += 1
                            sted = val
                            val += 1
                            return (g-grains_left)

        if steady_state == 0:
            total_array[0][val] = heights[0][0]
            val += 1
            avalanche_array_tran.append(relaxation)
        else:
            if grains_left > 0:
                total_array[0][val] = heights[0][0]
                avalanche_array_steady.append(relaxation)
                val += 1

def main():

    height , avalanches, tran_heights, slope = oslo_algorithm(5, 40)
    print(avalanches)

if __name__ == "__main__":
    main()