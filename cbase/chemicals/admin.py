from django.contrib import admin
from .models import *


admin.site.register(Premesis, DepartmentRole, Department, UserRole, User, Shelf, ChemicalBox,
                    Supplier, Chemical, Bottle)
