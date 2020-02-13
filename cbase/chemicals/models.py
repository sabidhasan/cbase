from django.db import models

class Premesis(models.Model):
  # A physical premesis representing a company site 
  # has_many departments
  name = models.CharField(max_length=50)
  address = models.CharField(max_length=250)
  building = models.CharField(max_length=100)


class DepartmentRole(models.Model):
  # Roles for departments
  # has_many departments
  name = models.CharField(max_length=50)


class Department(models.Model):
  # A department, eg. Chemistry, Biology, Animal Studies
  # has_many shelves and users
  name = models.CharField(max_length=50)
  floor = models.IntegerField(default=1)
  contact_information = models.CharField(max_length=250, null=True)
  role = models.ForeignKey(DepartmentRole, on_delete=models.PROTECT)
  premesis = models.ForeignKey(Premesis, on_delete=models.PROTECT)


class UserRole(models.Model):
  # Roles for users: head, chemist, biologist, safety, shipping/recieving
  # has_many users
  name = models.CharField(max_length=50)


class User(models.Model):
  # A user of chemicals
  # has_many chemical_boxes
  name = models.CharField(max_length=50)
  email = models.EmailField(null=True)
  phone = models.CharField(max_length=25, null=True)
  department = models.ForeignKey(Department, on_delete=models.CASCADE)
  role = models.ForeignKey(UserRole, on_delete=models.CASCADE)


class Shelf(models.Model):
  # A department consists of many shelves, which contain chemical_boxes, which have chemicals
  # has_many chemical_boxes
  name = models.CharField(max_length=50)
  description = models.CharField(max_length=200, null=True)
  department = models.ForeignKey(Department, on_delete=models.PROTECT)


class ChemicalBox(models.Model):
  # A chemical_box contains many chemicals
  # has_many bottles
  name = models.CharField(max_length=50)
  shelf = models.ForeignKey(Shelf, on_delete=models.PROTECT)


class Supplier(models.Model):
  # A supplier represents a vendor for a bottle of a chemical
  # has_many bottles
  name = models.CharField(max_length=25)
  website = models.URLField()


class Chemical(models.Model):
  # Represents a chemical entity
  # has_many bottles
  name = models.CharField(max_length=100)
  cas = models.CharField(max_length=25, null=True)
  synonyms = models.TextField(null=True)
  # Hill notation formula
  formula = models.CharField(max_length=250, null=True)
  molecular_weight = models.FloatField(default=0.0)
  density = models.FloatField(default=0.0)
  smiles = models.CharField(max_length=500, null=True)
  description = models.CharField(max_length=250, null=True)


class Bottle(models.Model):
  # A physical bottle of a chemical
  barcode = models.CharField(max_length=50)
  size = models.FloatField()
  size_unit = models.CharField(max_length=10, choices=[('g', 'grams'), ('ml', 'milliliters')])
  amount_remaining = models.FloatField()
  cost = models.IntegerField(default=0)
  # Whether the bottle is fresh and in a usable state
  usable = models.BooleanField(default=True)
  supplier = models.ForeignKey(Supplier, on_delete=models.PROTECT)
  box = models.ForeignKey(ChemicalBox, on_delete=models.PROTECT)
  chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
