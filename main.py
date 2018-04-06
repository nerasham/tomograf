from skimage import data, io, color
import numpy as np
import skimage
import math
from numpy import pi
from dicom.dataset import Dataset, FileDataset
import datetime
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import tkinter as tk
from tkinter import messagebox
import string
from tkinter import ttk
from tkinter import filedialog
import matplotlib.pyplot as plt
import dicom
from scipy import linalg


def bresenham_line(x1, y1, x2, y2, wi, hi):
    line = []
    x = x1
    y = y1
    if x1 < x2:
        xi = 1
        dx = x2 - x1;
    else:
        xi = -1
        dx = x1 - x2

    if y1 < y2:
        yi = 1
        dy = y2 - y1
    else:
        yi = -1
        dy = y1 - y2
    if x > -wi and x < wi and y > -hi and y < hi:
        line.append((x, y))

    if dx > dy:
        ai = (dy - dx) * 2
        bi = dy * 2
        d = bi - dx
        while x != x2:
            if d >= 0:
                x += xi
                y += yi
                d += ai
            else:
                d += bi
                x += xi
            if x > -wi and x < wi and y > -hi and y < hi:
                line.append((x, y))
    else:
        ai = (dx - dy) * 2
        bi = dx * 2
        d = bi - dy
        while y != y2:
            if d >= 0:
                x += xi
                y += yi
                d += ai
            else:
                d += bi
                y += yi
            if x > -wi and x < wi and y > -hi and y < hi:
                line.append((x, y))
    return line
alliter=[]

def createSinogram(picture, n_steps, l, n_detectors):
    # promień koła opisanego
    r = math.sqrt(picture.shape[0] ** 2 + picture.shape[1] ** 2)
    print(picture.shape[0])
    # środek koła opisanego
    w = picture.shape[0] // 2
    h = picture.shape[0] // 2
    # wynikowa lista
    sinogram=[]
    for step in range(n_steps):

        alpha = 2 * pi / n_steps * step
        # współrzędne emitera na kole
        emiter_point = np.array([r * math.cos(alpha), r * math.sin(alpha)]).astype(int)
        # współrzedne emitera na obrazku
        rays = []

        for i in range(n_detectors):

            # współrzędne detektora na kole
            detector_point = np.array([r * math.cos(alpha + pi + l / 2 - i * (l / (n_detectors - 1))), r * math.sin(alpha + pi + l / 2 - i * (l / (n_detectors - 1)))]).astype(int)
            # addytywne pochłanianie
            rays2 = bresenham_line(emiter_point[0], emiter_point[1], detector_point[0], detector_point[1], w, h)
            ray = 0
            if len(rays2) > 0:
                  for point in rays2:
                      ray += picture[point[0] + w - 1, point[1] + h - 1]
                  ray /= len(rays2)
            rays.append(ray)

        sinogram.append(rays)

    return sinogram


def createImage(sin, n_steps, l, n_detectors, image_width):
    # promień koła opisanego
    r = math.sqrt(image_width ** 2 + image_width ** 2)
    # środek koła opisanego
    w = image_width // 2
    h = image_width // 2
    # wynikowa lista
    image = np.zeros((image_width, image_width))
    print(len(sin))
    for step in range(n_steps):
        alpha = 2 * pi / n_steps * step
        # współrzędne emitera na kole
        emiter_point = np.array([r * math.cos(alpha), r * math.sin(alpha)]).astype(int)
        # współrzedne emitera na obrazku

        for i in range(n_detectors):

            # współrzędne detektora na kole
            detector_point = np.array([r * math.cos(alpha + pi + l / 2 - i * (l / (n_detectors - 1))),
                                       r * math.sin(alpha + pi + l / 2 - i * (l / (n_detectors - 1)))]).astype(int)

            # addytywne pochłanianie
            for point in bresenham_line(emiter_point[0], emiter_point[1], detector_point[0], detector_point[1], w, h):
                """print(point)"""
                image[point[0] + w - 1, point[1] + h - 1] += sin[step, i]

    return image


def kernel():
    length = 15
    middle = length // 2
    kern = np.zeros(length)
    for i in range(len(kern)):
        if middle == i:
            kern[i] = 1
        elif i % 2 != middle % 2:
            kern[i] = -4. / ((pi * abs(i - middle)) ** 2)
    return kern


def filter_1D(line, kernel):
    result = []

    middle = len(kernel) // 2
    for i in range(len(line)):
        sum = 0.
        tmp = 0
        if i < middle:
            temp = middle - i
            while temp < len(kernel):
                # print(temp)
                sum += kernel[temp] * line[tmp]
                tmp += 1
                temp += 1
            result.append(sum)
        elif i >= middle and i < len(line) - middle:
            while tmp < len(kernel):
                sum += kernel[tmp] * line[i + tmp - middle]
                tmp += 1
            result.append(sum)
        else:
            temp = middle + (len(line) - i)
            while tmp < temp:
                sum += kernel[tmp] * line[i + tmp - temp]
                tmp += 1
            result.append(sum)
    return result


def sinogram_filter(sin):
    k = kernel()
    result = []
    for i in range(len(sin)):
        temp = filter_1D(sin[i], k)
        result.append(temp)
    wyn = np.array(result)

    print()
    return wyn

def LoadImages(images,image,n_steps, n_detectors, l):
    images.changePlot((2,3,1),image)
    images.image_width=image.shape[0]
    images.n_detectors=n_detectors
    images.l=l
    images.per_steps=n_steps//20
    sinogram=createSinogram(image,n_steps,l,n_detectors)
    images.sin=sinogram
    sinogram=np.array(sinogram)
    sinogram/=np.max(sinogram)
    print(sinogram)
    images.changePlot((2, 3, 2), sinogram)
    im1=createImage(sinogram,n_steps,l,n_detectors,image.shape[0])
    im1/=np.max(im1)
    images.changePlot((2, 3, 3), im1)
    filtered_sin=sinogram_filter(sinogram)
    im2=createImage(filtered_sin,n_steps,l,n_detectors,image.shape[0])
    im2/=np.max(im2)
    images.changePlot((2, 3, 5), im2)
    for i in range(len(filtered_sin)):
        for j in range(len(filtered_sin[i])):
            if filtered_sin[i][j]!=0:
                filtered_sin[i][j]=abs(filtered_sin[i][j])
    filtered_sin/=np.max(filtered_sin)
    images.changePlot((2, 3, 4), filtered_sin)
    images.canvas.show();

def readDicomData(filename):
    file = readFileDICOM(filename)

    name = str(file.PatientName)
    name = name.split('^')
    lastName = name[0]
    firstName = name[1]
    id = str(file.PatientID)
    birthday = str(file.PatientBirthDate)
    birthday = birthday[0:4] + "-" + birthday[4:6] + "-" + birthday[6:8]
    gender = str(file.PatientSex)
    date = str(file.StudyDate)  # Study date
    time = str(file.StudyTime)   # Study time
    date = date[0:4]+"-"+date[4:6]+"-"+date[6:8]
    time = time[0:2]+":"+time[2:4]+":"+time[4:6]
    image = file.pixel_array
    return lastName, firstName, id,birthday, gender, date, time, image


def readFileDICOM(filename):
    return dicom.read_file(filename)

def saveFileDICOM(filename, patientName, patientId, sex, birth, imageArray, transpose=False):
    meta = Dataset()
    SOPClassUID = "1.2.840.10008.5.1.4.1.1.2" # sop class UID dla obrazow CT
    meta.MediaStorageSOPClassUID = SOPClassUID  # Wygenerowany unikalny UID
    date=datetime.datetime.now().strftime('%Y%m%d') # Obecny czas
    time=datetime.datetime.now().strftime('%H%M%S.%f') # Obecny czas
    randId = SOPClassUID + "."+date+time # Wygenerowany unikalny UID
    meta.MediaStorageSOPInstanceUID = randId # Wygenerowany unikalny UID
    meta.ImplementationClassUID = randId+"."+"1" # Wygenerowany unikalny UID

    allData = FileDataset(filename, {}, file_meta=meta, preamble=b"\0"*128) # Utworzenie obiektu DICOM
    allData.PatientName = patientName   # Imie pacjenta
    allData.PatientID=patientId # Id pacjenta
    allData.PatientsBirthDate = birth # Data urodzenia pacjenta
    allData.PatientSex = sex # Plec pacjenta
    allData.is_little_endian=True
    allData.is_implicit_VR=True
    allData.ContentDate = date  # Czas utworzenia pliku (YYYY:MM:DD)
    allData.StudyDate = date    # Czas ostatniego otworzenia obrazu (YYYY-MM-DD)
    allData.StudyTime = time    # Czas ostatniego otworzenia obrazu (HH:MM:SS)
    allData.ContentTime=time    # Czas utworzenia pliku (HH:MM:SS)
    allData.StudyInstanceUID = randId+"."+"2"   # Wygenerowany unikalny UID
    allData.SeriesInstanceUID = randId+"."+"3"   # Wygenerowany unikalny UID
    allData.SOPInstanceUID = randId+"."+"4"   # Wygenerowany unikalny UID
    allData.SOPClassUID = "CT."+date+time   # Wygenerowany unikalny UID

    allData.SamplesPerPixel = 1 # Liczba kanałów. 1 - dla skali szarosci
    allData.PhotometricInterpretation = "MONOCHROME2" # MONOCHROE - obraz jest w skali szarości, 2 - maksymalna wartosc wskazuje kolor bialy
    allData.PixelRepresentation = 0 # 0 - wartosci sa tylko dodatnie (unsigned) 1 - wartosci sa tez ujemne
    allData.HighBit = 15    # Najważniejszy bit w pliku z obrazem
    allData.BitsStored = 16 # Liczba bitow na jedna wartosc w obrazie
    allData.BitsAllocated = 16  # Liczba bitow na jedna wartosc ktora jest zaalokowana dla obrazu
    allData.SmallestImagePixelValue = b'\\x00\\x00' # Wskazanie minimalnej wartosci dla kanalu
    allData.LargestImagePixelValue = b'\\xff\\xff'  # Wskazanie maksymalnej wartosci dla kanalu
    allData.Rows = imageArray.shape[1]  # Liczba wierszy
    allData.Columns = imageArray.shape[0]   # Liczba kolumn
    if imageArray.dtype != np.uint16:   # Sprawdzenie czy wartosci sa w przedziale [0,255]
        imageArray = skimage.img_as_uint(imageArray)    # Zamiana na wartosci w przedziale [0,255]
        if transpose == True:   # Zamiana wierszy i kolumn (opcjonalne)
            allData.Rows = imageArray.shape[0]
            allData.Columns = imageArray.shape[1]
    allData.PixelData = imageArray.tostring()   # Zapisanie obrazu
    allData.save_as(filename) # Zapisanie pliku na dysku

def ret_dif_bet_2_inputs(input1,input2):
    if (len(input1) != len(input2)) or (len(input1[0]) != len(input2[0])): raise NameError(
        "Rozmiary tablicy nie są równe")
    val = 0
    for x in range(0, len(input1)):
        for y in range(0, len(input1[0])):
            val += (input1[x][y] - input2[x][y])** 2
    return np.sqrt(val)

def savetest(x,data,filename,labelname):
    plik=open(filename+'.txt','w')
    plik.write(str(x)+"\n\n")
    plik.write(str(data)+"\n\n")
    plik.close()
    plt.gcf().clear()
    plt.plot(x, data, '--bo')
    plt.ylabel("wariancja")
    plt.xlabel(labelname)
    plt.savefig(filename + ".pdf")

def choosedate(x,y,z,data,filename,labelname):
    if(len(x)==1 and len(y)==1 and len(z)>1):
        savetest(z,data[0,0,:],filename,labelname)
    if(len(x)==1 and len(y)>1 and len(z)==1):
        savetest(y,data[0,:,0],filename,labelname)
    elif(len(x)>0 and len(y)==0 and len(z)==0):
        savetest(x,data[:,0,0],filename,labelname)

def test_iterations(n_steps, l, n_detectors,filename,labelname,filtr):
    image=data.imread("os.png",as_grey=True)
    image/=np.max(image)
    per_steps=np.arange(0,361,20)
    per_steps[len(per_steps)-1]=360
    result=np.zeros(len(per_steps))
    for  i,j in enumerate(per_steps):
        if j==0:
            j+=1
        sin=createSinogram(image,j,l,n_detectors)
        sin2=sin.copy()
        sin2/=np.max(sin2)
        if filtr==True:
            sin2=sinogram_filter(sin2)
            invImage=createImage(sin2,j,l,n_detectors,200)
        else:
            invImage = createImage(sin2, j, l, n_detectors, 200)

        cpy=invImage.copy()
        if filtr==True:
            for x in range(len(cpy)):
                for y in range(len(cpy[x])):
                    if cpy[x][y] != 0:
                        cpy[x][y] = abs(cpy[x][y])
        cpy/=np.max(cpy)
        result[i]=ret_dif_bet_2_inputs(image,cpy)

        savetest(per_steps, result, filename, labelname)

l_param=np.arange(0,361,36)
l_param[0]=1
l_param[10]=359
l_param=l_param*(pi/180)

n_det = np.arange(50,401,35)
n_st=np.arange(0,361,36)
n_st[10]=359
n_st[0]=1

def testsAlgorithm(n_steps, nr_l, n_detectors,filename,labelname,filtr):
    image = data.imread("os.png", as_grey=True)
    image /= np.max(image)
    result=np.zeros((len(nr_l),len(n_detectors),len(n_steps)))
    invImage=None
    for l,lval in enumerate(nr_l):
        for d,dval in enumerate(n_detectors):
            for s,sval in enumerate(n_steps):
                print(lval,dval,sval)
                sin=createSinogram(image,sval,lval,dval)
                sin = np.array(sin)
                sin/=np.max(sin)
                if filtr==True:
                    sin2=sinogram_filter(sin)
                    sin2 = np.array(sin2)
                    sin2/=np.max(sin2)
                    invImage=createImage(sin2,sval,lval,dval,200)
                    invImage/=np.max(invImage)
                else:
                    invImage=createImage(sin,sval,lval,dval,200)
                    invImage /= np.max(invImage)
                result[l,d,s]=ret_dif_bet_2_inputs(image,invImage)
    choosedate(nr_l, n_detectors,n_steps, result, filename, labelname)

class Graphical(tk.Tk):
    val=0
    def __init__(self):

        tk.Tk.__init__(self)
        tk.Tk.wm_title(self, "CT")
        self.path = "No file"
        """print(l_param)
        print(n_det)
        print(n_st)"""
        self.n_steps = ["N Krokow", 1, 45]
        self.n_detectors = ["N Detektorow", 25, 400]
        self.prudence_angle=["Rozwartosc [rad]",30,180]
        control, leftpannel,middlepannel,rightpannel = self.Menu()
        self.images = AllImages()
        container = tk.Frame(self)
        container.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        container.tkraise()
        canvas = FigureCanvasTkAgg(self.images, container)
        self.images.canvas = canvas
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

    def Menu(self):
        padx = 20
        pady = 5
        textWidth = 20
        value = 360
        globall = 0.
        global_n_steps=0
        control = tk.Frame(self)
        control.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        leftpannel = tk.Frame(control)
        leftpannel.grid(row=1, column=1, sticky=tk.W)
        n_steps = ttk.Labelframe(leftpannel, text=self.n_steps[0])
        n_steps.grid(row=1, column=1, padx=padx, pady=pady, sticky=tk.W)
        slider3 = tk.Scale(n_steps, from_=20, to=360, orient=tk.HORIZONTAL, resolution=20, length=140,command=lambda value: self.val2(int(value)))
        slider3.pack(padx=padx, pady=pady)
        slider3.set(360)
        n_detectors = ttk.Labelframe(leftpannel, text=self.n_detectors[0])
        n_detectors.grid(row=2, column=1, padx=padx, pady=pady, sticky=tk.W)
        self.n_detectorsText = tk.Text(n_detectors, width=textWidth, height=n_detectors.winfo_height())
        self.n_detectorsText.pack(padx=padx, pady=pady)
        self.n_detectorsText.insert("1.0", "200")
        prudence_angle = ttk.Labelframe(leftpannel, text=self.prudence_angle[0])
        prudence_angle.grid(row=3,column=1,padx=padx,pady=pady)
        slider2 = tk.Scale(prudence_angle, from_=0., to=2*pi,orient=tk.HORIZONTAL,resolution = 0.02,length=140,command=lambda value: self.val1(float(value)))
        slider2.pack(padx=padx, pady=pady)
        slider2.set(pi/2)

        middlepannel = tk.Frame(control)
        middlepannel.grid(row=1, column=2)
        self.DICOMframe=ttk.Labelframe(middlepannel, text="DICOM")
        self.DICOMframe.grid(row=1,column=2,padx=padx,pady=pady,sticky=tk.W)
        lastName = ttk.Labelframe(self.DICOMframe, text="Nazwisko")
        lastName.grid(row=1, column=1, padx=padx, pady=pady, sticky=tk.W)
        self.lastNameDicom = tk.Text(lastName, width=textWidth, height=lastName.winfo_height())
        self.lastNameDicom.pack()
        Name = ttk.Labelframe(self.DICOMframe, text="Imie")
        Name.grid(row=2, column=1, padx=padx, pady=pady, sticky=tk.W)
        self.NameDicom = tk.Text(Name, width=textWidth, height=Name.winfo_height())
        self.NameDicom.pack()
        sex = ttk.Labelframe(self.DICOMframe, text="Plec")
        sex.grid(row=3, column=1, padx=padx, pady=pady, sticky=tk.W)
        self.sexx = tk.StringVar()
        self.sexx = ttk.Combobox(sex, justify=tk.CENTER, state="readonly", width=5,height=sex.winfo_height(), textvariable=self.sexx)
        self.sexx.pack()
        self.sexx["values"] = ('M', 'K')
        self.sexx.current(0)
        birth = ttk.Labelframe(self.DICOMframe, text="Data ur.")
        birth.grid(row=1, column=2, padx=padx, pady=pady, sticky=tk.W)
        self.birthh = tk.Text(birth, width=textWidth, height=birth.winfo_height())
        self.birthh.pack()
        idDicomFrame = ttk.Labelframe(self.DICOMframe, text="Id")
        idDicomFrame.grid(row=2, column=2, padx=padx, pady=pady, sticky=tk.W)
        self.idDicom = tk.Text(idDicomFrame, width=textWidth, height=idDicomFrame.winfo_height())
        self.idDicom.pack()
        dateDicomFrame = ttk.Labelframe(self.DICOMframe, text="Data(YYYY-MM-DD)")
        dateDicomFrame.grid(row=3, column=2, padx=padx, pady=pady, sticky=tk.W)
        self.dateDicom = tk.Label(dateDicomFrame, text="Data")
        self.dateDicom.pack()
        timeDicomFrame = ttk.Labelframe(self.DICOMframe, text="Czas(HH:MM:SS)")
        timeDicomFrame.grid(row=1, column=3, padx=padx, pady=pady, sticky=tk.W)
        self.timeDicom = tk.Label(timeDicomFrame, text="Czas")
        self.timeDicom.pack()
        Savebutton = ttk.Button(self.DICOMframe, text="Zapisz plik", command=lambda: self.saveFile())
        Savebutton.grid(row=3, column=3, padx=padx, pady=pady)

        rightpannel = tk.Frame(control)
        rightpannel.grid(row=1, column=3)
        fileFrame = ttk.Labelframe(rightpannel, text="Opcje ")
        fileFrame.grid(row=1, column=3, padx=padx, pady=pady, sticky=tk.W)
        self.filterVar = tk.BooleanVar()
        Loadbutton = ttk.Button(fileFrame, text="Otwórz plik", command=lambda: self.loadFile())
        Loadbutton.grid(row=1, column=1, padx=padx, pady=pady)
        self.LoadLabel = ttk.Label(fileFrame, text=self.path)
        self.LoadLabel.grid(row=1, column=2, padx=padx, pady=pady, sticky=tk.W)

        GenerateButton = ttk.Button(fileFrame, text="Odtworz obraz", command=lambda: self.generate())
        GenerateButton.grid(row=2, column=1, padx=padx, pady=pady)
        testbutton = ttk.Button(fileFrame, text="Test statystyczny", command=lambda: self.testgraph())
        testbutton.grid(row=3, column=1, padx=padx, pady=pady)

        sliderLabel = ttk.Labelframe(rightpannel, text="Zmiana obrazu")
        sliderLabel.grid(row=4, column=3,padx=padx, pady=pady, sticky=tk.W)

        slider = tk.Scale(sliderLabel, from_=1, to=20,orient=tk.HORIZONTAL,length=170,command=lambda value: self.images.swapImage(int(value)-1))
        slider.set(20)
        slider.pack(padx=padx, pady=pady)

        return control, leftpannel,middlepannel,rightpannel
    def testgraph(self):
        test_iterations(360,pi/2,200,'test1','n krokow',False)
        test_iterations(360,pi/2,200,'test2','n krokow',True)

        testsAlgorithm([180],[pi/2],n_det,'test3','n detektorow',False)
        testsAlgorithm([180],[pi/2],n_det,'test4','n detektorow',True)
        testsAlgorithm(n_st, [pi / 2], [200],'test5','n krokow', False)
        testsAlgorithm(n_st, [pi / 2], [200],'test6','n krokow', True)
        testsAlgorithm([180], l_param, [200], 'test7','rozwartosc katowa',False)
        testsAlgorithm([180], l_param, [200],'test8','rozwartosc katowa', True)

        """print("krok 5")
        result9 = testsAlgorithm([180],l_param, n_det, False)
        plik9 = open('test9.txt', 'w')
        plik9.write(str(l_param) + "\n\n")
        plik9.write(str(n_det) + "\n\n")
        plik9.write(str(result9) + "\n\n")
        plik9.close()
        result10 = testsAlgorithm([180],l_param, n_det, True)
        plik10 = open('test10.txt', 'w')
        plik10.write(str(l_param) + "\n\n")
        plik10.write(str(n_det) + "\n\n")
        plik10.write(str(result10) + "\n\n")
        plik10.close()"""

    def generate(self):
        detectors = self.n_detectorsText.get("1.0", "end-1c")
        LoadImages(self.images,self.images.plots[(2,3,1)][1],n_steps=self.global_n_steps,n_detectors=int(detectors),l=self.globall)

    def val1(self,value):
        self.globall=value

    def val2(self,value):
        self.global_n_steps=value

    def saveFile(self):
        lastName = self.lastNameDicom.get("1.0", "end-1c")
        if lastName == "":
            messagebox.showinfo(title="Nazwisko", message=("Nazwisko nie zostalo podane"))
            return
        Name = self.NameDicom.get("1.0", "end-1c")
        if Name == "":
            messagebox.showinfo(title="Imie", message=("Imie nie zostalo podane"))
            return
        sex = self.sexx.get()
        birth = self.birthh.get("1.0", "end-1c")
        if len(birth) != 8:
            messagebox.showinfo(title="Data Urodzenia", message=("Zla data urodzenia"))
            return
        id = self.idDicom.get("1.0", "end-1c")
        if id == "":
            messagebox.showinfo(title="Pusty id", message=("Id pacjenta jest pusty"))
            return
        filePath = filedialog.asksaveasfilename(filetypes=(("DICOM files", "*.dcm"),
                                                           ("All files", "*.*")))
        if filePath == '':
            return
        filePath = filePath.split(".")[0]
        extension = ".dcm"
        saveFileDICOM(filePath + "Original" + extension, lastName + '^' + Name, id, sex, birth,
                           self.images.plots[(2, 3, 1)][1])
        saveFileDICOM(filePath + "Sinogram" + extension, lastName + '^' + Name, id, sex, birth,
                           self.images.plots[(2, 3, 2)][1], transpose=True)

        saveFileDICOM(filePath + "Imagewithsin" + extension, lastName + '^' + Name, id, sex, birth,
                  self.images.plots[(2, 3, 3)][1], transpose=True)
        saveFileDICOM(filePath + "Sinogramwithfiltered" + extension, lastName + '^' + Name, id, sex, birth,
                      self.images.plots[(2, 3, 4)][1], transpose=True)
        saveFileDICOM(filePath + "Imagewithfiltered" + extension, lastName + '^' + Name, id, sex, birth,
                      self.images.plots[(2, 3, 5)][1], transpose=True)

    def loadFile(self):
        self.path = filedialog.askopenfilename()
        if self.path == () or self.path == "":
            self.path = ""
            return
        self.LoadLabel.config(text=self.path)
        fileExtension = self.path.split(".")
        if len(fileExtension) < 2 or not fileExtension[1] == "dcm":
            newImage = data.imread(self.path, as_grey=True)
            self.images.changePlot((2,3, 1), newImage)
            return

        lastNamez, firstNamez, id,birthdayz, sexz, datez, timez, newImage = readDicomData(self.path)
        self.images.changePlot((2, 3, 1), newImage)
        self.lastNameDicom.delete("1.0", tk.END)
        self.lastNameDicom.insert("1.0", lastNamez)
        self.NameDicom.delete("1.0", tk.END)
        self.NameDicom.insert("1.0", firstNamez)
        self.idDicom.delete("1.0", tk.END)
        self.idDicom.insert("1.0", id)
        self.dateDicom.config(text=datez)
        self.timeDicom.config(text=timez)
        self.birthh.delete("1.0", tk.END)
        self.birthh.insert("1.0", birthdayz)
        sex = sexz.upper()
        if sex == "M":
            self.sexx.current(0)
        elif sex == "F":
            self.sexx.current(1)


class AllImages(Figure):
    plots={}
    canvas=None
    image_width=0
    n_detectors=0
    l=0
    per_steps=0
    sin=[]
    def __init__(self, figsize=(10, 10), dpi=100):
        Figure.__init__(self, figsize=figsize, dpi=dpi)
        titles = [{"title": "Orginalny obraz"}, {"title": "sinogram"},
                  {"title": "IRadon"},{"title": "Sinogram po filtrze"},{"title": "IRadon po filtrze"}]
        for i in range(1,6):
            plt=self.add_subplot(2,3,i)
            if "title" in titles[i-1]:
                 plt.set_title(titles[i-1]["title"])
            plt.set_xticks([])
            plt.set_yticks([])
            self.plots[(2,3,i)]=[plt,None]

    def changePlot(self, plot, image, cmap="gray"):
        self.plots[plot][1] = image
        self.plots[plot][0].imshow(image,cmap=cmap)
        self.canvas.show()

    def swapImage(self, number):
        if number < 1 or number > 22:
            return
        if self.sin!=[]:
            sy=np.asarray(self.sin[:number*self.per_steps])
            sy/=np.max(sy)
            self.changePlot((2,3,2),sy)
            print(self.l)
            print(self.n_detectors)
            print(self.image_width)
            imag=createImage(sy,number*self.per_steps,self.l,self.n_detectors,self.image_width)
            imag/=np.max(imag)
            self.changePlot((2, 3, 3),imag)
            fs = sinogram_filter(sy)
            im3 = createImage(fs, number*self.per_steps,self.l, self.n_detectors, self.image_width)
            im3 /= np.max(im3)
            self.changePlot((2, 3, 5), im3)
            for i in range(len(fs)):
                for j in range(len(fs[i])):
                    if fs[i][j] != 0:
                        fs[i][j] = abs(fs[i][j])
            fs /= np.max(fs)
            self.changePlot((2, 3, 4), fs)


def main():

    app =Graphical()
    app.mainloop()


if __name__ == "__main__":
    main()