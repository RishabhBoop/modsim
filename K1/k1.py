# ------------------------------------------------
# Aufgabe K1
# ------------------------------------------------

# ------------------------------------------------
import numpy as np
import scipy
import matplotlib.pyplot as plt
# ------------------------------------------------

def aufgabe1():
    print("---------- Aufgabe 1 ----------")
    a = 3
    # b = [2, 3]
    b = np.array([2, 3])
    c = np.arange(1, 79, 2)

    # variablen ausgeben
    print("Variable erzeugt: a =", a)
    print("Variable erzeugt: a =", b)
    print("Variable erzeugt: a =", c)

    # variablen in eine datei schreiben
    np.savez("K1/test.npz", a=a, b=b, c=c)

    print("Variablen in Datei 'K1/test.npz' gespeichert.")

    # variablen löschen
    del a
    del b
    del c
    print("Variablen gelöscht.")

    # verifizieren dass die variablen gelöscht wurden
    print("Verifizierung dass die Variablen wirklich gelöscht wurden:")

    try:
        print(a)
    except NameError:
        print("Variable 'a' ist nicht definiert")

    try:
        print(b)
    except NameError:
        print("Variable 'b' ist nicht definiert")
    try:
        print(c)
    except NameError:
        print("Variable 'c' ist nicht definiert")


    # variablen aus der datei wieder einlesen
    data = np.load("K1/test.npz")
    print(data)
    a = data["a"]
    b = data["b"]
    c = data["c"]

    # verifizieren dass die variablen wieder eingelesen wurden
    print("Wieder eingelesene Variable a =", a)
    print("Wieder eingelesene Variable b =", b)
    print("Wieder eingelesene Variable c =", c)

    print("---------- Ende Aufgabe 1 ----------")

    return

# ------------------------------------------------

def aufgabe2():
    print("---------- Aufgabe 2 ----------")

    magic = np.loadtxt("K1/magic_matrix.csv", delimiter=",")
    print("Form der Matrix:", magic.shape)

    print("Summe aller Elemente:", np.sum(magic))
    print("Maximum aller Elemente:", np.max(magic))
    p_max = np.unravel_index(np.argmax(magic), magic.shape)
    print("Position des Maximums aller Elemente:", np.argmax(magic), "bei Position [", p_max[0], ",", p_max[1], "]")
    p_min = np.unravel_index(np.argmin(magic), magic.shape)
    print("Position des Minimums aller Elemente:", np.argmin(magic), "bei Position [", p_min[0], ",", p_min[1], "]")

    v1 = magic[:, 2].reshape(magic.shape[0], 1)
    print(v1)

    print("Summe aller Elemente die größer als 40 sind:", np.sum(magic[magic > 40]))

    pos_1 = np.ravel_multi_index((3, 4), magic.shape)
    print("pos 3,4:", pos_1)
    pos_2 = np.ravel_multi_index((7, 4), magic.shape)
    print("pos 7,4:", pos_2)
    v2 = magic.ravel()[pos_1 : pos_2 + 1]
    print(v2)
    print("---------- Ende Aufgabe 2 ----------")
    return


def aufgabe3():
    print("---------- Aufgabe 3 ----------")
    print("...")
    print("---------- Ende Aufgabe 3 ----------")
    return

def aufgabe4():
    print("---------- Aufgabe 4 ----------")
    # load .mat file
    matfile = scipy.io.loadmat("K1/mkl.mat")
    # print(matfile)

    # extract values
    xWerte = matfile.get("x_werte").squeeze()
    yWerte = matfile.get("y_werte").squeeze()
    xWerte = np.array(xWerte)
    yWerte = np.array(yWerte)
    # print(xWerte)
    # print(yWerte)

    # interpolate values
    interpol = scipy.interpolate.interp1d(xWerte, yWerte)
    x_new = np.arange(min(xWerte), max(xWerte), 0.1)
    y_new = interpol(x_new)

    # plot matfile content
    plt.plot(xWerte, yWerte, color='blue')
    plt.plot(x_new, y_new, color='black', linestyle='--')
    print("Showing plot...")
    plt.show()
    print("Plot closed.")
    print("---------- Ende Aufgabe 4 ----------")
    return

aufgabe_choice = input("Welche Aufgabe möchten Sie ausführen? ")
if aufgabe_choice == "1":
    aufgabe1()
elif aufgabe_choice == "2":
    aufgabe2()
elif aufgabe_choice == "3":
    aufgabe3()
elif aufgabe_choice == "4":
    aufgabe4()
else:
    print("Ungültige Eingabe. Bitte gültige Zahl eingeben.")
