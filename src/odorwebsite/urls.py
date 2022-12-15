"""ServerOdors URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.views.defaults import server_error

from . import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.urls import path, include
#from rest_framework import routers

from odor.views import index, ketcher, OdorWebSite, OdorWebSite_Home, OdorWebSite_About, OdorWebSite_Search, \
    OdorWebSite_Predict, OdorWebSite_Contact, OdorWebSite_OlfactoryReceptor_template, OdorWebSite_Chemical_template, \
    OdorWebSite_search_chem_or_odor, OdorWebSite_Odor_template, OdorWebSite_design, test, OdorWebSite_docking_chem_or
#router = routers.DefaultRouter()
#router.register(r'todos', views.TodoView, 'todo')

urlpatterns = [
    #path('admin/', admin.site.urls),
    #path('api/', include(router.urls)),
    #path('guillaume/', server_error, name="er"),
    path('index', index, name="index"),
    path('', OdorWebSite_Home, name="OdorWebSite_Home"),
    path('about', OdorWebSite_About, name="OdorWebSite_About"),
    path('search-structure', OdorWebSite_Search, name="OdorWebSite_Search"),
    path('predict', OdorWebSite_Predict, name="OdorWebSite_Predict"),
    path('contact', OdorWebSite_Contact, name="OdorWebSite_Contact"),
    path('or/<idOlfactoryReceptors>', OdorWebSite_OlfactoryReceptor_template, name="OdorWebSite_OR"),
    path('chemical/<chem_id>', OdorWebSite_Chemical_template, name="OdorWebSite_Chemical"),
    path('odor/<odor_id>', OdorWebSite_Odor_template, name="OdorWebSite_Odor"),
    path('search', OdorWebSite_search_chem_or_odor, name="search_chem_or_odor"),
    path('docking', OdorWebSite_docking_chem_or, name="docking_chem_or"),
    path('ketcher', ketcher, name="ketcher"),
    path('design', OdorWebSite_design, name="design"),
    #path("", include("Server.urls"))
]
# TemplateView.as_view(template_name="scripts/ketcher-2.4.0/index.html"),
#if settings.DEBUG:
#    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)