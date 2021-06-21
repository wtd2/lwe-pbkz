import os

mark = [0, 60, 65, 70, 80, 85, 90, 100, 110, 120, 120, 130, 140, 150, 160, 170, 180, 190, 200]

os.system("rm verify ; make verify")
score = 0
for i in range(1, 19):
    if os.system("./verify %d > vrf/%d.txt 2>> vrf/%d.txt" % (i, i, i)) == 0:
        print("%d\t\033[1;32mPASS\033[0m\t%d" % (i, mark[i]))
        score += mark[i]
    else:
        print("%d\t\033[1;31mFAIL\033[0m\t%d" % (i, mark[i]))
print("score: %d/%d" % (score, sum(mark)))