require('ateles_lineareuler')
name = 'ateles_lineareuler'


input = {
         --read = './restart/simulation_lastHeader.lua',
         read = 'restart/ateles_euler_lastHeader.lua',

         subsampling = {
                         levels = math.floor(math.log(degree+1)/math.log(2))+1,
                         --levels = 1,
                         projection = 'QLegendrePoint',
                       },

  }


output = {
            --folder = './harvester/',
            folder = 'harvester/',

           {


            format = 'VTU',

            binary = true,
            requestedData = {
                             variable = {
                                          { name = 'density', },
                                          { name = 'momentum',},
                                          { name = 'energy',  },
                                          { name = 'pressure', ncomponents = 1 },
                                          { name = 'velocity', ncomponents = 3 },
                                        },
                            },
            vrtx = {
                   },
           }

         }
