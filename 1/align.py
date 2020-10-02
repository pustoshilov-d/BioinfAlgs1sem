import numpy as np


class NidlmanVunsh:
    def __init__(self, str1, str2, points):
        self.InDel = points["InDel"]
        self.Swap = points["Swap"]
        self.Match = points["Match"]

        self.str1 = str1
        self.str2 = str2

        self.xLen = len(str1)
        self.yLen = len(str2)

        self.matrix = np.array([1000] * (self.xLen + 1) * (self.yLen + 1)).reshape(
            self.xLen + 1, self.yLen + 1
        )

        self.result = self.culcMatrix(self.xLen, self.yLen)
        self.score = self.result["score"]
        self.align = self.result["align"]

    def culcMatrix(self, i, j):

        curAlign = ["", ""]

        if i == 0:
            if j == 0:
                curScore = 0
                curAlign = ["", ""]

            else:
                previous = self.culcMatrix(i, j - 1)
                curScore = previous["score"] + self.InDel
                curAlign[0] = previous["align"][0] + "-"
                curAlign[1] = previous["align"][1] + self.str2[j - 1]

        elif j == 0:
            previous = self.culcMatrix(i - 1, j)
            curScore = previous["score"] + self.InDel
            curAlign[0] = previous["align"][0] + self.str1[i - 1]
            curAlign[1] = previous["align"][1] + "-"

        else:
            if self.str1[i - 1] == self.str2[j - 1]:
                previous = self.culcMatrix(i - 1, j - 1)
                curScore = previous["score"] + self.Match
                curAlign[0] = previous["align"][0] + self.str1[i - 1]
                curAlign[1] = previous["align"][1] + self.str2[j - 1]

            else:
                previous = {
                    "swap": self.culcMatrix(i - 1, j - 1),
                    "in": self.culcMatrix(i, j - 1),
                    "del": self.culcMatrix(i - 1, j),
                }

                previous["swap"]["score"] += self.Swap
                previous["in"]["score"] += self.InDel
                previous["del"]["score"] += self.InDel

                curScore = max(
                    previous["swap"]["score"],
                    previous["in"]["score"],
                    previous["del"]["score"],
                )

                if previous["swap"]["score"] == curScore:
                    curAlign[0] = previous["swap"]["align"][0] + self.str1[i - 1]
                    curAlign[1] = previous["swap"]["align"][1] + self.str2[j - 1]
                elif previous["in"]["score"] == curScore:
                    curAlign[0] = previous["in"]["align"][0] + "-"
                    curAlign[1] = previous["in"]["align"][1] + self.str2[j - 1]
                elif previous["del"]["score"] == curScore:
                    curAlign[0] = previous["del"]["align"][0] + self.str1[i - 1]
                    curAlign[1] = previous["del"]["align"][1] + "-"

        self.matrix[i, j] = curScore

        # print(i,j)
        # print(self.matrix)
        # print(curScore)
        # print(curAlign)

        return {"score": curScore, "align": curAlign}


if __name__ == "__main__":

    string1 = "магия"
    string2 = "химия"
    points = {"InDel": -1, "Swap": -1, "Match": 0}

    alg = NidlmanVunsh(string1, string2, points)
    print("Matrix:")
    print(alg.matrix)
    print("Score: ", alg.score)
    print("Alignment:")
    print(alg.align[0])
    print(alg.align[1])
    print()

    string1 = "GFAMYW"
    string2 = "ALYFAM"
    points = {"InDel": -2, "Swap": -0, "Match": 3}

    alg = NidlmanVunsh(string1, string2, points)
    print("Matrix:")
    print(alg.matrix)
    print("Score: ", alg.score)
    print("Alignment:")
    print(alg.align[0])
    print(alg.align[1])
    print()

    string1 = "GFAMYW"
    string2 = "ALYFAM"
    points = {"InDel": -0, "Swap": -1, "Match": 1}
    alg = NidlmanVunsh(string1, string2, points)
    print("Matrix:")
    print(alg.matrix)
    print("Score: ", alg.score)
    print("Alignment:")
    print(alg.align[0])
    print(alg.align[1])
    print()
