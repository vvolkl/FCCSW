
#include "DetCommon/DetUtils.h"

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/IActsExtension.hpp"
#include "Acts/Plugins/Digitization/CartesianSegmentation.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "DD4hep/DetFactoryHelper.h"

using dd4hep::Volume;
using dd4hep::DetElement;
using dd4hep::xml::Dimension;
using dd4hep::PlacedVolume;

namespace det {

std::shared_ptr<const Acts::DigitizationModule>
rectangleDigiModuleXZ(double halflengthX, double halflengthZ, double thickness, double gridSizeX, double gridSizeZ) {
  // convert to ACTS units
  double scalor = Acts::units::_cm;
  halflengthX *= scalor;
  halflengthZ *= scalor;
  thickness *= scalor;
  auto bounds = std::make_shared<const Acts::RectangleBounds>(halflengthX, halflengthZ);
  // the Acts segmentation of the DigitizationModule
  size_t bins0 = (gridSizeX != 0) ? (2 * halflengthX) / (gridSizeX * scalor) : 0;
  size_t bins1 = (gridSizeZ != 0) ? (2 * halflengthZ) / (gridSizeZ * scalor) : 0;

  std::shared_ptr<const Acts::CartesianSegmentation> actsSegmentation =
      std::make_shared<const Acts::CartesianSegmentation>(bounds, bins0, bins1);

  // finally create the digitization module
  // @todo set lorentz angle
  return (std::make_shared<const Acts::DigitizationModule>(actsSegmentation, thickness, 1, 0));
};

static dd4hep::Ref_t createTkLayoutTrackerBarrel(dd4hep::Detector& lcdd,
                                                 dd4hep::xml::Handle_t xmlElement,
                                                 dd4hep::SensitiveDetector sensDet) {
  // shorthands
  dd4hep::xml::DetElement xmlDet = static_cast<dd4hep::xml::DetElement>(xmlElement);
  Dimension dimensions(xmlDet.dimensions());
  // get sensitive detector type from xml
  dd4hep::xml::Dimension sdTyp = xmlElement.child(_Unicode(sensitive));
  // sensitive detector used for all sensitive parts of this detector
  sensDet.setType(sdTyp.typeStr());

  // definition of top volume
  // has min/max dimensions of tracker for visualization etc.
  std::string detectorName = xmlDet.nameStr();
  DetElement topDetElement(detectorName, xmlDet.id());
  Acts::ActsExtension::Config barrelConfig;
  barrelConfig.isBarrel = true;
  // detElement owns extension
  Acts::ActsExtension* detWorldExt = new Acts::ActsExtension(barrelConfig);
  topDetElement.addExtension<Acts::IActsExtension>(detWorldExt);
  dd4hep::Tube topVolumeShape(dimensions.rmin(), dimensions.rmax(), (dimensions.zmax() - dimensions.zmin()) * 0.5);
  Volume topVolume(detectorName, topVolumeShape, lcdd.air());
  topVolume.setVisAttributes(lcdd.invisible());

  // counts all layers - incremented in the inner loop over repeat - tags
  unsigned int layerCounter = 0;
  double integratedModuleComponentThickness = 0;
  double phi = 0;
  // loop over 'layer' nodes in xml
  dd4hep::xml::Component xLayers = xmlElement.child(_Unicode(layers));
  for (dd4hep::xml::Collection_t xLayerColl(xLayers, _U(layer)); nullptr != xLayerColl; ++xLayerColl) {
    std::vector<dd4hep::xml::Component> rodLists;
    dd4hep::xml::Component xLayer = static_cast<dd4hep::xml::Component>(xLayerColl);
    dd4hep::xml::Component xRods = xLayer.child("rods");
    rodLists.push_back(xRods.child("rodOdd"));
    rodLists.push_back(xRods.child("rodEven"));
    rodLists.push_back(xRods.child("rodOddTilted", false));
    rodLists.push_back(xRods.child("rodEvenTilted", false));
    dd4hep::xml::Component xModulePropertiesOdd = rodLists[0].child("moduleProperties");
    dd4hep::Tube layerShape(xLayer.rmin(), xLayer.rmax(), dimensions.zmax());
    Volume layerVolume("layer", layerShape, lcdd.material("Air"));
    // layerVolume.setVisAttributes(lcdd.invisible());
    PlacedVolume placedLayerVolume = topVolume.placeVolume(layerVolume);
    placedLayerVolume.addPhysVolID("layer", layerCounter);
    DetElement lay_det(topDetElement, "layer" + std::to_string(layerCounter), layerCounter);
    Acts::ActsExtension::Config layConfig;
    layConfig.isLayer = true;
    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    layConfig.axes = "XzY";  // correct translation of local x axis in dd4hep to local x axis in acts
    // detElement owns extension
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension(layConfig);
    lay_det.addExtension<Acts::IActsExtension>(layerExtension);
    lay_det.setPlacement(placedLayerVolume);
    dd4hep::xml::Component xModuleComponentsOdd = xModulePropertiesOdd.child("components");
    integratedModuleComponentThickness = 0;
    int moduleCounter = 0;
    Volume moduleVolume;
    // store the materials of each tracking module for tracking geometry
    std::vector<std::pair<dd4hep::Material, double>> compMaterials;
    // go through all components to collect the material to later attach it to the modoule
    for (dd4hep::xml::Collection_t xModuleComponentOddColl(xModuleComponentsOdd, _U(component));
         nullptr != xModuleComponentOddColl;
         ++xModuleComponentOddColl) {
      dd4hep::xml::Component xModuleComponentOdd = static_cast<dd4hep::xml::Component>(xModuleComponentOddColl);
      compMaterials.push_back(
          std::make_pair(lcdd.material(xModuleComponentOdd.materialStr()), xModuleComponentOdd.thickness()));
    }

    for (dd4hep::xml::Collection_t xModuleComponentOddColl(xModuleComponentsOdd, _U(component));
         nullptr != xModuleComponentOddColl;
         ++xModuleComponentOddColl) {
      dd4hep::xml::Component xModuleComponentOdd = static_cast<dd4hep::xml::Component>(xModuleComponentOddColl);
      double moduleWidth = 0.5 * xModulePropertiesOdd.attr<double>("modWidth");
      double moduleThickness = 0.5 * xModuleComponentOdd.thickness();
      double moduleLength = 0.5 * xModulePropertiesOdd.attr<double>("modLength");

      moduleVolume = Volume("module",
                            dd4hep::Box(moduleWidth, moduleThickness, moduleLength),
                            lcdd.material(xModuleComponentOdd.materialStr()));

      // Create digitization module
      // with readout given by layer
      auto digiModule = rectangleDigiModuleXZ(moduleWidth, moduleLength, moduleThickness, 0.025, 0.05);

      moduleVolume.setVisAttributes(lcdd.invisible());
      unsigned int nPhi = xRods.repeat();
      for (unsigned int phiIndex = 0; phiIndex < nPhi; phiIndex += 2) {
        double lX = 0;
        double lY = 0;
        double lZ = 0;
        for (auto currentComp : rodLists) {
          if (currentComp.ptr()) {

            phi = 2 * M_PI * static_cast<double>(phiIndex) / static_cast<double>(nPhi);
            for (dd4hep::xml::Collection_t xModuleColl(currentComp.child(_Unicode(modules)), _U(module));
                 nullptr != xModuleColl;
                 ++xModuleColl) {
              dd4hep::xml::Component xModule = static_cast<dd4hep::xml::Component>(xModuleColl);
              double currentPhi = atan2(xModule.Y(), xModule.X());
              double componentOffset = integratedModuleComponentThickness -
                  0.5 * xModulePropertiesOdd.attr<double>("modThickness") + 0.5 * xModuleComponentOdd.thickness();
              dd4hep::Translation3D offsetOnly(cos(currentPhi) * componentOffset, sin(currentPhi) * componentOffset, 0);
              lX = xModule.X();  // + cos(currentPhi) * componentOffset;
              lY = xModule.Y();  // + sin(currentPhi) * componentOffset;
              lZ = xModule.Z();
              dd4hep::Translation3D moduleOffset(lX, lY, lZ);
              dd4hep::Transform3D lTrafo(dd4hep::RotationZ(atan2(lY, lX) + 0.5 * M_PI), moduleOffset);
              dd4hep::RotationZ lRotation(phi);
              double thetaTilt = xModule.attr<double>("thetaTilt");
              dd4hep::RotationX lRotation_thetaTilt(-1 * thetaTilt);
              PlacedVolume placedModuleVolume =
                  layerVolume.placeVolume(moduleVolume, lRotation * lTrafo * lRotation_thetaTilt * offsetOnly);
              if (xModuleComponentOdd.isSensitive()) {
                placedModuleVolume.addPhysVolID("module", moduleCounter);
                moduleVolume.setSensitiveDetector(sensDet);
                DetElement mod_det(lay_det, "module" + std::to_string(moduleCounter), moduleCounter);

                // add extension to hand over material
                Acts::ActsExtension* moduleExtension = new Acts::ActsExtension(compMaterials, digiModule);
                mod_det.addExtension<Acts::IActsExtension>(moduleExtension);

                mod_det.setPlacement(placedModuleVolume);
                ++moduleCounter;
              }
            }
          }
        }
      }
      integratedModuleComponentThickness += xModuleComponentOdd.thickness();
    }
    ++layerCounter;
  }
  Volume motherVol = lcdd.pickMotherVolume(topDetElement);
  PlacedVolume placedGenericTrackerBarrel = motherVol.placeVolume(topVolume);
  placedGenericTrackerBarrel.addPhysVolID("system", topDetElement.id());
  topDetElement.setPlacement(placedGenericTrackerBarrel);
  return topDetElement;
}
}  // namespace det

DECLARE_DETELEMENT(TkLayoutBrlTracker, det::createTkLayoutTrackerBarrel)
