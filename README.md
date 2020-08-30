# Ten-Bar Truss Homework

這是最佳化實驗室暑假作業，作業內容為SOLab_manual_v1.3.pdf中的2.2.1 Ten-bar Truss的部分。

這次作業主要要解ten-bar的問題，再將其以matlab求出最佳解，最後用latex撰寫一份pdf繳交。
推導過程與最後答案撰寫於Ten_Bar_Homework.Pdf。
最佳化的matlab程式碼為tenbar.m。

Latex於overleaf撰寫，其檔案至於Latex資料夾中。

## About Matlab Code
第一部我先設定fmincon的參數，接著用fmincon解ten-bar的問題，再以最佳解帶入公式計算出題目所要求的答案。
### function objfun
這是最佳化的objective function，在本題目為所有桿件的總重量
### function nonlcon
這是最佳化的非線性boundary condition撰寫的函數，輸出c為線性不等式，ceq為線性等式。其中線性不等式的撰寫方式是該函式須小於零。
本題邊界條件有11個，包含了10個應力與1個位移條件。
### function answer
fmincon求出的r1 r2最佳解於程式中代碼為x，此函數由最佳解計算出題幹所要求的位移disp、應力stress、反作用力force。
