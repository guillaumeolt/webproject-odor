# Generated by Django 4.1.2 on 2022-10-13 09:23

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('odor', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='chemicalsodors',
            name='SMILE',
            field=models.CharField(max_length=250, unique=True),
        ),
    ]
