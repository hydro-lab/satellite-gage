for (i in (2:f)){
    if (alng_per[i,2]>threepeak){
        if (alng_per[i-1,2]<threepeak){
            i1 <- alng_per[i-1,1]
            i2 <- alng_per[i,1]
            j1 <- alng_per[i-1,2]
            j2 <- alng_per[i,2]
            n <- threepeak    
            RDB <- ((n-(j1))*((i2-i1)/(j2-j1))+i1)
            break
            }
    }
}
for (i in 1:(f-1)){
    if (alng_per[f-i,2]>threepeak){        #expressing the index such that when i = 1, f, and when i = 2, f-1.
        if (alng_per[f-i+1,2]<threepeak){
            i1 <- alng_per[f-i+1,1]
            i2 <- alng_per[f-1,1]
            j1 <- alng_per[f-1+1,2]
            j2 <- alng_per[f-1,2]
            n <- threepeak    
            LDB <- ((n-(j1))*((i2-i1)/(j2-j1))+i1)
            break
        }
    }
}

width <- LDB-RDB #gives width in meters
