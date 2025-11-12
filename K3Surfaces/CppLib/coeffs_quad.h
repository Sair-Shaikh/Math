#define ABC \
    A = mult[mult[mult[mult[y_0][y_0]][y_0]][y_0]][y_1] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_0]][y_2]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_2]][y_2]][y_2]][y_2] ^ \
   mult[mult[mult[mult[y_1][y_1]][y_1]][y_2]][y_2] ^ \
   mult[mult[mult[mult[y_1][y_1]][y_2]][y_2]][y_2] ^ \
   mult[mult[mult[mult[y_2][y_2]][y_2]][y_2]][y_2], \
    B = mult[mult[mult[mult[y_0][y_0]][y_0]][y_0]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_0]][y_1]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_2]][y_2]][y_2], \
    C = mult[mult[mult[mult[y_0][y_0]][y_0]][y_0]][y_0] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_0]][y_0]][y_1] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_0]][y_0]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_1]][y_2]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_1]][y_1]][y_2]][y_2] ^ \
   mult[mult[mult[mult[y_1][y_1]][y_1]][y_1]][y_1] ^ \
   mult[mult[mult[mult[y_1][y_1]][y_1]][y_1]][y_2] ^ \
   mult[mult[mult[mult[y_1][y_1]][y_2]][y_2]][y_2]


#define ABC2 \
    A = mult[mult[y_0][y_2]][y_2] ^ \
   mult[mult[y_1][y_1]][y_1] ^ \
   mult[mult[y_2][y_2]][y_2], \
    B = mult[mult[y_0][y_1]][y_2] ^ \
   mult[mult[y_2][y_2]][y_2], \
    C = mult[mult[y_0][y_2]][y_2] ^ \
   mult[mult[y_1][y_1]][y_1] ^ \
   mult[mult[y_1][y_1]][y_2] ^ \
   mult[mult[y_2][y_2]][y_2]