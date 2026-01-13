EXPORT_PLOT = True
EXPORT_ALL = True


def export_chapter1():
    import K1.k1 as k1

    k1.EXPORT_PLOT = EXPORT_PLOT
    k1.aufgabe4()
    print("Kapitel 1 wurde exportiert.")


def export_chapter2():
    import K2.k2 as k2

    k2.EXPORT_PLOT = EXPORT_PLOT
    if EXPORT_ALL:
        k2.aufgabe1()
        k2.aufgabe2()
        k2.aufgabe3()
    k2.aufgabe5()
    print("Kapitel 2 wurde exportiert.")


def export_chapter3():
    import K3.k3 as k3

    k3.EXPORT_PLOT = EXPORT_PLOT
    if EXPORT_ALL:
        k3.aufgabe1()
        k3.aufgabe2()
        k3.aufgabe3()
    k3.aufgabe5()
    print("Kapitel 3 wurde exportiert.")


def export_chapter4():
    import K4.k4 as k4

    k4.EXPORT_PLOT = EXPORT_PLOT
    if EXPORT_ALL:
        k4.aufgabe1_2()
        k4.aufgabe3()
        k4.aufgabe4()
    k4.aufgabe5()
    k4.aufgabe6()
    print("Kapitel 4 wurde exportiert.")


def export_chapter5():
    import K5.k5 as k5

    k5.EXPORT_PLOT = EXPORT_PLOT
    if EXPORT_ALL:
        k5.aufgabe1()
        k5.aufgabe2()
        k5.aufgabe3()
    k5.aufgabe4()
    print("Kapitel 5 wurde exportiert.")


def export_chapter6():
    import K6.k6 as k6

    k6.EXPORT_PLOT = EXPORT_PLOT
    if EXPORT_ALL:
        k6.aufgabe1(export_filename="k6_a1.png")
        k6.aufgabe2()
    k6.aufgabe4()
    print("Kapitel 6 wurde exportiert.")


print("Hallo, dieses Skript exportiert die plots der verschiedenen Kapitel.")

export_all_input = input("Möchten Sie alle Kapitel exportieren? (y/n): ").strip().lower()
if export_all_input == "y":
    EXPORT_ALL = True
elif export_all_input == "n":
    EXPORT_ALL = False
else:
    print("Ungültige Eingabe, verwende Standardwert (alle Aufgaben)")
    EXPORT_ALL = True

chapters = input(
    "Welche Kapitel möchten Sie exportieren?\nGeben Sie die Kapitelnnummern durch Kommas getrennt ein, ENTER für alle (z.B. 1,3,5): "
)
if chapters.strip() == "":
    chapters_list = ["1", "2", "3", "4", "5", "6"]
else:   
    chapters_list = [ch.strip() for ch in chapters.split(",")]

for chapter in chapters_list:
    if chapter == "1":
        export_chapter1()
    elif chapter == "2":
        export_chapter2()
    elif chapter == "3":
        export_chapter3()
    elif chapter == "4":
        export_chapter4()
    elif chapter == "5":
        export_chapter5()
    elif chapter == "6":
        export_chapter6()
    else:
        print(f"Kapitel {chapter} ist nicht verfügbar zum Exportieren.")
