import numpy as np
import gmsh
import vtk
import math
import os
import random


# Класс расчётной сетки
class CalcMesh:

    # Конструктор сетки, полученной из stl-файла
    def __init__(self, nodes_coords, tetrs_points):
        # 3D-сетка из расчётных точек
        # Пройдём по узлам в модели gmsh и заберём из них координаты
        self.nodes = np.array([nodes_coords[0::3],nodes_coords[1::3],nodes_coords[2::3]])

        # Модельная скалярная величина распределена как-то вот так
        self.smth = np.power(self.nodes[0, :], 2) + np.power(self.nodes[1, :], 2)

        # Тут может быть скорость, но сейчас здесь нули
        self.velocity = np.zeros(shape=(3, int(len(nodes_coords) / 3)), dtype=np.double)
        self.sizedata=int(len(nodesCoord) / 3)
        # Пройдём по элементам в модели gmsh
        self.tetrs = np.array([tetrs_points[0::4],tetrs_points[1::4],tetrs_points[2::4],tetrs_points[3::4]])
        self.tetrs -= 1

    # Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    def move(self, tau, vect):
        # По сути метод просто двигает все точки c их текущими скоростями
        #self.velocity = vect
        self.nodes += (vect * tau)

    def rotate(self, vect):
        self.nodes += vect
        self.velocity=vect


    # Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    def snapshot(self, snap_number):
        # Сетка в терминах VTK
        unstructuredGrid = vtk.vtkUnstructuredGrid()
        # Точки сетки в терминах VTK
        points = vtk.vtkPoints()

        # Скалярное поле на точках сетки
        smth = vtk.vtkDoubleArray()
        smth.SetName("smth")

        # Векторное поле на точках сетки
        vel = vtk.vtkDoubleArray()
        vel.SetNumberOfComponents(3)
        vel.SetName("vel")

        # Обходим все точки нашей расчётной сетки
        # Делаем это максимально неэффективным, зато наглядным образом
        for i in range(0, len(self.nodes[0])):
            # Вставляем новую точку в сетку VTK-снапшота
            points.InsertNextPoint(self.nodes[0,i], self.nodes[1,i], self.nodes[2,i])
            # Добавляем значение скалярного поля в этой точке
            smth.InsertNextValue(self.smth[i])
            # Добавляем значение векторного поля в этой точке
            vel.InsertNextTuple((self.velocity[0,i], self.velocity[1,i], self.velocity[2,i]))

        # Грузим точки в сетку
        unstructuredGrid.SetPoints(points)

        # Присоединяем векторное и скалярное поля к точкам
        unstructuredGrid.GetPointData().AddArray(smth)
        unstructuredGrid.GetPointData().AddArray(vel)

        # А теперь пишем, как наши точки объединены в тетраэдры
        # Делаем это максимально неэффективным, зато наглядным образом
        for i in range(0, len(self.tetrs[0])):
            tetr = vtk.vtkTetra()
            for j in range(0, 4):
                tetr.GetPointIds().SetId(j, self.tetrs[j,i])
            unstructuredGrid.InsertNextCell(tetr.GetCellType(), tetr.GetPointIds())

        # Создаём снапшот в файле с заданным именем
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputDataObject(unstructuredGrid)
        writer.SetFileName("apovorot" + str(snap_number) + ".vtu")
        writer.Write()


# Теперь придётся немного упороться:
# (а) построением сетки средствами gmsh,
# (б) извлечением данных этой сетки в свой код.
gmsh.initialize()

# Считаем STL
try:
    path = os.path.dirname(os.path.abspath(__file__))
    gmsh.merge(os.path.join(path, 'grill.stl'))
except:
    print("Could not load STL mesh: bye!")
    gmsh.finalize()
    exit(-1)

# Восстановим геометрию
angle = 40
forceParametrizablePatches = False
includeBoundary = True
curveAngle = 180
gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary, forceParametrizablePatches, curveAngle * math.pi / 180.)
gmsh.model.mesh.createGeometry()

# Зададим объём по считанной поверхности
s = gmsh.model.getEntities(2)
l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])
gmsh.model.geo.addVolume([l])

gmsh.model.geo.synchronize()

# Зададим мелкость желаемой сетки
f = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(f, "F", "5")
gmsh.model.mesh.field.setAsBackgroundMesh(f)

# Построим сетку
gmsh.model.mesh.generate(3)

# Теперь извлечём из gmsh данные об узлах сетки
nodeTags, nodesCoord, parametricCoord = gmsh.model.mesh.getNodes()

# И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
GMSH_TETR_CODE = 4
tetrsNodesTags = None
elementTypes, elementTags, elementNodeTags = gmsh.model.mesh.getElements()
for i in range(0, len(elementTypes)):
    if elementTypes[i] != GMSH_TETR_CODE:
        continue
    tetrsNodesTags = elementNodeTags[i]

if tetrsNodesTags is None:
    print("Can not find tetra data. Exiting.")
    gmsh.finalize()
    exit(-2)

print("The model has %d nodes and %d tetrs" % (len(nodeTags), len(tetrsNodesTags) / 4))

# На всякий случай проверим, что номера узлов идут подряд и без пробелов
for i in range(0, len(nodeTags)):
    # Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
    assert (i == nodeTags[i] - 1)
# И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
assert(len(tetrsNodesTags) % 4 == 0)

# TODO: неплохо бы полноценно данные сетки проверять, да
###################################################################################
def calc (n, angle): #ось - прямая, массив из 3 чисел [n1, n2, n3] и вращение относительно о
    v1 = []
    v2 = []
    v3 = []
    o=[0,1,0]
    for i in range(int(len(nodesCoord) / 3)):
        err2 = np.zeros(3)
        
        r = [nodes[0][i]-21, nodes[1][i], nodes[2][i]-1.5]

        modul_n0 = math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])

        n[0] = n[0]/modul_n0
        n[1] = n[1]/modul_n0
        n[2] = n[2]/modul_n0
        scal = r[0] * n[0] + r[1] * n[1] + r[2] * n[2]
        modul_n = math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])

        err = [0, 0, 0]
        err[0] = r[0] - (scal*n[0])/modul_n
        err[1] = r[1] - (scal * n[1]) / modul_n
        err[2] = r[2] - (scal * n[2]) / modul_n

        modul_rad = math.sqrt(err[0]*err[0] + err[1]*err[1] + err[2]*err[2])
        err1 = np.array(err)/modul_rad

        if err != [0, 0, 0]:
            s = np.cross(np.array(err), np.array(n)/modul_n)
            modul_s = math.sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2])
            err2 = s / modul_s
            err3 = (-math.sin(angle/2)*err1 + math.cos(angle/2)*err2)* modul_rad*math.sin(angle/2)*2
        else:
            err3 = np.zeros(3)
        #print(max(nodes[2]))
        v1.append(err3[0])
        v2.append(err3[1])
        v3.append(err3[2])


    return(np.array([v1, v2, v3]))
###################################################################################################################################################


mesh = CalcMesh(nodesCoord, tetrsNodesTags)
nodes = mesh.nodes
tau = 1
#v = np.zeros(shape=(3, int(len(nodesCoord) / 3)), dtype=np.double)


mesh.snapshot(0)
for i in range(1, int(150)):
    nodes = mesh.nodes
    print("nodes1:", i, "-", nodes[0][0], nodes[1][0], nodes[2][0])
    #print(math.sqrt(((nodes[1][0]-nodes[1][1])*(nodes[1][0]-nodes[1][1])) + ((nodes[0][0]-nodes[0][1])*(nodes[0][0]-nodes[0][1])) + ((nodes[2][0]-nodes[2][1])*(nodes[2][0]-nodes[2][1]))))
    mesh.rotate(calc([0, 1, 0], math.pi/6))
    mesh.snapshot(i)



gmsh.finalize()
