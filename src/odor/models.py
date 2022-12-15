from django.db import models
from django.db import connection
# Create your models here.


class Todo(models.Model):
    title = models.CharField(max_length=120)
    description = models.TextField()
    completed = models.BooleanField(default=False)

    def _str_(self):
        return self.title


class ChemicalsOdors(models.Model):
    idChemicals = models.IntegerField(primary_key=True)
    Name = models.CharField(max_length=450)
    IUPAC_name = models.CharField(max_length=450)
    SMILE = models.CharField(max_length=250, unique=True)
    Pubchem_CID = models.CharField(max_length=450)
    Mixture = models.SmallIntegerField()
    Synonyms = models.CharField(max_length=800)
    CAS = models.CharField(max_length=45)
    Molecular_Weight = models.CharField(max_length=45)
    Molecular_Formula = models.CharField(max_length=45)
    InChi = models.CharField(max_length=400, default='None')
    InChi_Key = models.CharField(max_length=45, default='None')
    class Meta:
        db_table = "Chemicals"

class OlfactoryReceptors(models.Model):
    idOlfactoryReceptors = models.IntegerField(primary_key=True)
    GeneName = models.CharField(max_length=45)
    idUniprot = models.CharField(max_length=45)
    Synonym = models.CharField(max_length=45)
    Species = models.CharField(max_length=45)
    Sequence = models.CharField(max_length=5000)
    URL_3D_Structure = models.CharField(max_length=45)
    FileName_3D_Structure = models.CharField(max_length=45)
    class Meta:
        db_table = "OlfactoryReceptors"

class Smell_Percepts(models.Model):
    idSmell_Percepts = models.IntegerField(primary_key=True)
    Odor = models.CharField(max_length=45)
    class Meta:
        db_table = "Smell_Percepts"

def dictfetchall(cursor):
    "Return all rows from a cursor as a dict"
    columns = [col[0] for col in cursor.description]
    return [
        dict(zip(columns, row))
        for row in cursor.fetchall()
    ]

def my_custom_sql():
    cursor = connection.cursor()
    sql_query = "select mydb.Chemicals.Name, "\
                        "mydb.Chemicals.IUPAC_name, "\
                        "mydb.Chemicals.CAS, "\
                        "mydb.Chemicals.idChemicals, "\
                        "mydb.Chemicals.Pubchem_CID, "\
                        "mydb.Chemicals.SMILE, "\
                        "GROUP_CONCAT(distinct Smell_Percepts.Odor SEPARATOR ';') AS smell, "\
                        "GROUP_CONCAT(distinct OlfactoryReceptors.idOlfactoryReceptors SEPARATOR ';') as idOlfactoryReceptors, "\
                        "GROUP_CONCAT(distinct OlfactoryReceptors.GeneName SEPARATOR ';') as OlfRecept "\
                 "from mydb.Chemicals "\
                        "left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals "\
                        "left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts "\
                        "left join mydb.OlfactoryReceptors_has_Chemicals on mydb.Chemicals.idChemicals = mydb.OlfactoryReceptors_has_Chemicals.Chemicals_idChemicals "\
                        "left join mydb.OlfactoryReceptors on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.OlfactoryReceptors_has_Chemicals.OlfactoryReceptors_idOlfactoryReceptors "\
                     "GROUP BY idChemicals, Pubchem_CID, CAS, Name, IUPAC_name, SMILE"

    cursor.execute(sql_query);
    dic_db = dictfetchall(cursor)
    return(dic_db)

def my_custom_sql_with_human_homologue():
    cursor = connection.cursor()
    sql_query = "select mydb.Chemicals.idChemicals, "\
                        "mydb.Chemicals.Pubchem_CID, "\
                        "mydb.Chemicals.CAS, "\
                        "mydb.Chemicals.Name, "\
                        "mydb.Chemicals.IUPAC_name, "\
                        "mydb.Chemicals.SMILE, "\
                        "GROUP_CONCAT(distinct Smell_Percepts.Odor SEPARATOR ';') AS smell, "\
                        "GROUP_CONCAT(distinct NULLIF(CONCAT(COALESCE(OlfactoryReceptors.idOlfactoryReceptors,''),';',COALESCE(OlfRecept.idOlfactoryReceptors,'')),';') SEPARATOR ';') as idOlfactoryReceptors, "\
                        "GROUP_CONCAT(distinct NULLIF(CONCAT(COALESCE(OlfactoryReceptors.GeneName,''),';',COALESCE(OlfRecept.GeneName,'')),';') SEPARATOR ';') as OlfRecept "\
                 "from mydb.Chemicals "\
                        "left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals "\
                        "left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts "\
                        "left join mydb.OlfactoryReceptors_has_Chemicals on mydb.Chemicals.idChemicals = mydb.OlfactoryReceptors_has_Chemicals.Chemicals_idChemicals "\
                        "left join mydb.OlfactoryReceptors on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.OlfactoryReceptors_has_Chemicals.OlfactoryReceptors_idOlfactoryReceptors " \
                        "left join  mydb.RecepteurHomologueHumain on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.RecepteurHomologueHumain.OlfactoryReceptors_idOlfactoryReceptors "\
                        "left join mydb.OlfactoryReceptors OlfRecept on OlfRecept.idOlfactoryReceptors = mydb.RecepteurHomologueHumain.OlfactoryReceptors_idOlfactoryReceptors1 "\
                     "GROUP BY idChemicals, Pubchem_CID, CAS, Name, IUPAC_name, SMILE"

    cursor.execute(sql_query);
    dic_db = dictfetchall(cursor)
    return(dic_db)

def my_custom_sql_chem_id(chem_id):
    cursor = connection.cursor()
    sql_query = "select mydb.Chemicals.idChemicals, "\
                        "mydb.Chemicals.Pubchem_CID, "\
                        "mydb.Chemicals.CAS, "\
                        "mydb.Chemicals.Name, "\
                        "mydb.Chemicals.IUPAC_name, "\
                        "mydb.Chemicals.SMILE, " \
                        "GROUP_CONCAT(distinct Smell_Percepts.idSmell_Percepts SEPARATOR ';') AS idSmell, " \
                        "GROUP_CONCAT(distinct Smell_Percepts.Odor SEPARATOR ';') AS smell, "\
                        "GROUP_CONCAT(distinct OlfactoryReceptors.GeneName SEPARATOR ';') as OlfRecept "\
                 "from mydb.Chemicals "\
                        "left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals "\
                        "left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts "\
                        "left join mydb.OlfactoryReceptors_has_Chemicals on mydb.Chemicals.idChemicals = mydb.OlfactoryReceptors_has_Chemicals.Chemicals_idChemicals "\
                        "left join mydb.OlfactoryReceptors on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.OlfactoryReceptors_has_Chemicals.OlfactoryReceptors_idOlfactoryReceptors " \
                "WHERE mydb.Chemicals.idChemicals = %s "\
                     "GROUP BY idChemicals, Pubchem_CID, CAS, Name, IUPAC_name, SMILE "


    cursor.execute(sql_query, [str(chem_id)]);
    dic_db = dictfetchall(cursor)
    return(dic_db)
def my_custom_sql_odor_id(odor_id):
    cursor = connection.cursor()
    sql_query = "select  GROUP_CONCAT(distinct Chemicals.idChemicals SEPARATOR ';') as idChemicals, "\
                        "GROUP_CONCAT(distinct Chemicals.Pubchem_CID SEPARATOR ';') as Pubchem_CID, "\
                        "GROUP_CONCAT(distinct Chemicals.CAS SEPARATOR ';') as CAS, "\
                        "GROUP_CONCAT(distinct Chemicals.Name SEPARATOR ';') as Name, "\
                        "GROUP_CONCAT(distinct Chemicals.IUPAC_name SEPARATOR ';') as IUPAC_name, "\
                        "GROUP_CONCAT(distinct Chemicals.SMILE SEPARATOR ';') as SMILE, " \
                        "Smell_Percepts.idSmell_Percepts, " \
                        "Smell_Percepts.Odor, "\
                        "GROUP_CONCAT(distinct OlfactoryReceptors.GeneName SEPARATOR ';') as OlfRecept "\
                 "from mydb.Chemicals "\
                        "left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals "\
                        "left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts "\
                        "left join mydb.OlfactoryReceptors_has_Chemicals on mydb.Chemicals.idChemicals = mydb.OlfactoryReceptors_has_Chemicals.Chemicals_idChemicals "\
                        "left join mydb.OlfactoryReceptors on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.OlfactoryReceptors_has_Chemicals.OlfactoryReceptors_idOlfactoryReceptors " \
                "WHERE mydb.Smell_Percepts.idSmell_Percepts = %s "\
                     "GROUP BY idSmell_Percepts, Odor"


    cursor.execute(sql_query, [str(odor_id)]);
    dic_db = dictfetchall(cursor)
    return(dic_db)

def my_custom_sql_chem_predict():
    cursor = connection.cursor()
    sql_query = "select mydb.Chemicals.idChemicals, "\
                        "mydb.Chemicals.Pubchem_CID, "\
                        "mydb.Chemicals.CAS, "\
                        "mydb.Chemicals.Name, "\
                        "mydb.Chemicals.IUPAC_name, "\
                        "mydb.Chemicals.SMILE, "\
                        "GROUP_CONCAT(distinct Smell_Percepts.Odor SEPARATOR ';') AS smell, "\
                        "GROUP_CONCAT(distinct OlfactoryReceptors.GeneName SEPARATOR ';') as OlfRecept, "\
                        "training "\
                 "from mydb.Chemicals "\
                        "left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals "\
                        "left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts "\
                        "left join mydb.OlfactoryReceptors_has_Chemicals on mydb.Chemicals.idChemicals = mydb.OlfactoryReceptors_has_Chemicals.Chemicals_idChemicals "\
                        "left join mydb.OlfactoryReceptors on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.OlfactoryReceptors_has_Chemicals.OlfactoryReceptors_idOlfactoryReceptors " \
                        "left join mydb.Chemicals_has_ModelPredict on mydb.Chemicals.idChemicals = mydb.Chemicals_has_ModelPredict.Chemicals_idChemicals "\
                     "GROUP BY idChemicals, Pubchem_CID, CAS, Name, IUPAC_name, SMILE, training "


    cursor.execute(sql_query);
    dic_db = dictfetchall(cursor)
    return(dic_db)

def my_custom_sql_chem_predict():
    cursor = connection.cursor()
    sql_query = "select mydb.Chemicals.idChemicals, "\
                        "mydb.Chemicals.Pubchem_CID, "\
                        "mydb.Chemicals.CAS, "\
                        "mydb.Chemicals.Name, "\
                        "mydb.Chemicals.IUPAC_name, "\
                        "mydb.Chemicals.SMILE, "\
                        "GROUP_CONCAT(distinct Smell_Percepts.Odor SEPARATOR ';') AS smell, "\
                        "GROUP_CONCAT(distinct OlfactoryReceptors.GeneName SEPARATOR ';') as OlfRecept, "\
                        "training "\
                 "from mydb.Chemicals "\
                        "left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals "\
                        "left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts "\
                        "left join mydb.OlfactoryReceptors_has_Chemicals on mydb.Chemicals.idChemicals = mydb.OlfactoryReceptors_has_Chemicals.Chemicals_idChemicals "\
                        "left join mydb.OlfactoryReceptors on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.OlfactoryReceptors_has_Chemicals.OlfactoryReceptors_idOlfactoryReceptors " \
                        "left join mydb.Chemicals_has_ModelPredict on mydb.Chemicals.idChemicals = mydb.Chemicals_has_ModelPredict.Chemicals_idChemicals "\
                     "GROUP BY idChemicals, Pubchem_CID, CAS, Name, IUPAC_name, SMILE, training "


    cursor.execute(sql_query);
    dic_db = dictfetchall(cursor)
    return(dic_db)

def my_custom_sql_chem_get_odor(chem_id):
    cursor = connection.cursor()
    sql_query = "select  Smell_Percepts.idSmell_Percepts, " \
                        "Smell_Percepts.Odor "\
                 "from mydb.Chemicals "\
                        "left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals "\
                        "left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts "\
                "WHERE mydb.Chemicals.idChemicals = %s "


    cursor.execute(sql_query, [str(chem_id)])
    dic_db = dictfetchall(cursor)

    dic_db_k = dict()
    for k in dic_db:
        dic_db_k[k["idSmell_Percepts"]] = k["Odor"]
    return([dic_db_k])

def my_custom_sql_odor_get_chem(odor_id):
    cursor = connection.cursor()
    sql_query = "select  Chemicals.idChemicals, " \
                         "Chemicals.CAS, " \
                         "Chemicals.Pubchem_CID, " \
                         "Chemicals.IUPAC_name, " \
                         "Chemicals.SMILE, " \
                         "Chemicals.Name " \
                "from mydb.Chemicals "\
                        "left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals "\
                        "left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts "\
                "WHERE Smell_Percepts.idSmell_Percepts = %s ;"


    cursor.execute(sql_query, [str(odor_id)])
    dic_db = dictfetchall(cursor)

    dic_db_k = dict()
    for k in dic_db:
        dic_db_k[k["idChemicals"]] = k
    return([dic_db_k])

def my_custom_sql_or_get_chem(or_id):
    cursor = connection.cursor()
    sql_query = "select  Chemicals.idChemicals, " \
                         "Chemicals.CAS, " \
                         "Chemicals.Pubchem_CID, " \
                         "Chemicals.IUPAC_name, " \
                         "Chemicals.SMILE, " \
                         "Chemicals.Name " \
                "from mydb.Chemicals "\
                        "left join mydb.OlfactoryReceptors_has_Chemicals on mydb.Chemicals.idChemicals = mydb.OlfactoryReceptors_has_Chemicals.Chemicals_idChemicals "\
                        "left join mydb.OlfactoryReceptors on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.OlfactoryReceptors_has_Chemicals.OlfactoryReceptors_idOlfactoryReceptors "\
                "WHERE OlfactoryReceptors.idOlfactoryReceptors = %s ;"


    cursor.execute(sql_query, [str(or_id)])
    dic_db = dictfetchall(cursor)

    dic_db_k = dict()
    for k in dic_db:
        dic_db_k[k["idChemicals"]] = k
    return([dic_db_k])


def my_custom_sql_chem_get_odor_dic(chem_id):
    cursor = connection.cursor()
    sql_query = "select  Smell_Percepts.idSmell_Percepts, " \
                        "Smell_Percepts.Odor, "\
                        "Chemicals_has_Smell_Percepts.Sources, "\
                        "Chemicals_has_Smell_Percepts.odor_tag "\
                 "from mydb.Chemicals "\
                        "left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals "\
                        "left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts "\
                "WHERE mydb.Chemicals.idChemicals = %s "


    cursor.execute(sql_query, [str(chem_id)])
    dic_db = dictfetchall(cursor)

    dic_db_k = dict()
    for k in dic_db:
        dic_db_k[k["idSmell_Percepts"]] = k
    return(dic_db_k)

def my_custom_sql_chem_get_or_dic(chem_id):
    cursor = connection.cursor()
    sql_query = "select  OlfactoryReceptors.GeneName, " \
                        "OlfactoryReceptors.idOlfactoryReceptors, "\
                        "OlfactoryReceptors.idUniprot, "\
                        "OlfactoryReceptors.Synonym, "\
                        "OlfactoryReceptors.Species, "\
                        "OlfactoryReceptors_has_Chemicals.Sources "\
                 "from mydb.Chemicals "\
                        "left join mydb.OlfactoryReceptors_has_Chemicals on mydb.Chemicals.idChemicals = mydb.OlfactoryReceptors_has_Chemicals.Chemicals_idChemicals "\
                        "left join mydb.OlfactoryReceptors on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.OlfactoryReceptors_has_Chemicals.OlfactoryReceptors_idOlfactoryReceptors "\
                "WHERE GeneName IS NOT NULL and idUniprot IS NOT NULL and idChemicals = %s "
    print(sql_query)
    cursor.execute(sql_query, [str(chem_id)])
    dic_db = dictfetchall(cursor)

    dic_db_k = dict()
    for k in dic_db:
        dic_db_k[k["GeneName"]] = k
    return(dic_db_k)

"""

select  mydb.Chemicals.idChemicals,
        OlfactoryReceptors.GeneName,
        OlfactoryReceptors.idUniprot,
        OlfactoryReceptors.Synonym,
        OlfactoryReceptors.Species 
from mydb.Chemicals 
        left join mydb.OlfactoryReceptors_has_Chemicals on mydb.Chemicals.idChemicals = mydb.OlfactoryReceptors_has_Chemicals.Chemicals_idChemicals 
        left join mydb.OlfactoryReceptors on mydb.OlfactoryReceptors.idOlfactoryReceptors = mydb.OlfactoryReceptors_has_Chemicals.OlfactoryReceptors_idOlfactoryReceptors 
WHERE GeneName IS NOT NULL and idUniprot IS NOT NULL and idChemicals = 50;
"""
"""
def my_custom_sql_or_get_chem(or_id):


    select Odor, Sources, odor_tag 
from mydb.Chemicals 
        left join mydb.Chemicals_has_Smell_Percepts on mydb.Chemicals.idChemicals = mydb.Chemicals_has_Smell_Percepts.Chemicals_idChemicals 
        left join mydb.Smell_Percepts on mydb.Smell_Percepts.idSmell_Percepts = mydb.Chemicals_has_Smell_Percepts.Smell_Percepts_idSmell_Percepts 
WHERE Smell_Percepts.Odor !="" and idChemicals = 0;
"""