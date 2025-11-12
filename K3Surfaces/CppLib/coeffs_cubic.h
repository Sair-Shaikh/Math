// Sample coefficients: Key: 5-484431
#define ABCD \
    A = mult[mult[y_0][y_1]][y_1], \
    B = mult[mult[y_0][y_0]][y_0] ^ \
   mult[mult[y_0][y_0]][y_1] ^ \
   mult[mult[y_0][y_0]][y_2] ^ \
   mult[mult[y_0][y_1]][y_2] ^ \
   mult[mult[y_1][y_1]][y_1], \
    C = mult[mult[y_0][y_0]][y_0] ^ \
   mult[mult[y_0][y_2]][y_2] ^ \
   mult[mult[y_1][y_1]][y_1] ^ \
   mult[mult[y_1][y_1]][y_2], \
    D = mult[mult[y_0][y_0]][y_0] ^ \
   mult[mult[y_0][y_1]][y_1] ^ \
   mult[mult[y_0][y_2]][y_2] ^ \
   mult[mult[y_1][y_1]][y_1] ^ \
   mult[mult[y_2][y_2]][y_2]


#define ABCD2 \
    A = mult[mult[y_0][y_1]][y_1], \
    B = mult[mult[y_0][y_0]][y_0] ^ \
   mult[mult[y_0][y_0]][y_1] ^ \
   mult[mult[y_0][y_0]][y_2] ^ \
   mult[mult[y_0][y_1]][y_2] ^ \
   mult[mult[y_1][y_1]][y_1], \
    C = mult[mult[y_0][y_0]][y_0] ^ \
   mult[mult[y_0][y_2]][y_2] ^ \
   mult[mult[y_1][y_1]][y_1] ^ \
   mult[mult[y_1][y_1]][y_2], \
    D = mult[mult[y_0][y_0]][y_0] ^ \
   mult[mult[y_0][y_1]][y_1] ^ \
   mult[mult[y_0][y_2]][y_2] ^ \
   mult[mult[y_1][y_1]][y_1] ^ \
   mult[mult[y_2][y_2]][y_2]