// Sample coefficients: Key: 7-887005
#define ABCD \
    A = mult[mult[mult[mult[mult[y_0][y_0]][y_0]][y_0]][y_0]][y_1] ^ \
   mult[mult[mult[mult[mult[y_0][y_0]][y_0]][y_0]][y_0]][y_2] ^ \
   mult[mult[mult[mult[mult[y_0][y_0]][y_0]][y_1]][y_2]][y_2] ^ \
   mult[mult[mult[mult[mult[y_0][y_0]][y_1]][y_1]][y_2]][y_2] ^ \
   mult[mult[mult[mult[mult[y_0][y_1]][y_1]][y_1]][y_1]][y_1] ^ \
   mult[mult[mult[mult[mult[y_0][y_1]][y_1]][y_1]][y_1]][y_2] ^ \
   mult[mult[mult[mult[mult[y_1][y_1]][y_1]][y_1]][y_1]][y_1], \
    B = mult[mult[mult[mult[y_0][y_0]][y_0]][y_0]][y_1] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_0]][y_0]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_0]][y_1]][y_1] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_0]][y_1]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_0]][y_2]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_1]][y_1]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_0]][y_2]][y_2]][y_2] ^ \
   mult[mult[mult[mult[y_0][y_1]][y_1]][y_1]][y_1] ^ \
   mult[mult[mult[mult[y_1][y_1]][y_1]][y_1]][y_2], \
    C = mult[mult[mult[y_0][y_0]][y_0]][y_0] ^ \
   mult[mult[mult[y_0][y_0]][y_1]][y_1] ^ \
   mult[mult[mult[y_0][y_0]][y_2]][y_2] ^ \
   mult[mult[mult[y_0][y_1]][y_2]][y_2] ^ \
   mult[mult[mult[y_0][y_2]][y_2]][y_2] ^ \
   mult[mult[mult[y_1][y_1]][y_2]][y_2], \
    D = mult[mult[y_0][y_0]][y_0] ^ \
   mult[mult[y_0][y_0]][y_2] ^ \
   mult[mult[y_0][y_2]][y_2] ^ \
   mult[mult[y_2][y_2]][y_2]


#define ABCD2 \
    A = 1, \
    B = 0, \
    C = mult[mult[mult[y_0][y_1]][y_1]][y_2], \
    D = mult[mult[mult[mult[mult[y_0][y_0]][y_1]][y_2]][y_2]][y_2] ^ \
   mult[mult[mult[mult[mult[y_0][y_0]][y_2]][y_2]][y_2]][y_2] ^ \
   mult[mult[mult[mult[mult[y_1][y_1]][y_1]][y_1]][y_1]][y_2] ^ \
   mult[mult[mult[mult[mult[y_1][y_1]][y_1]][y_1]][y_2]][y_2] ^ \
   mult[mult[mult[mult[mult[y_1][y_1]][y_2]][y_2]][y_2]][y_2] ^ \
   mult[mult[mult[mult[mult[y_2][y_2]][y_2]][y_2]][y_2]][y_2]
