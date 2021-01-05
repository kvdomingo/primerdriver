import React from "react";
import { Image as CloudImage } from "cloudinary-react";
import { Link } from "react-router-dom";
import { PropTypes } from "prop-types";
import { Grid, Card, Row, Button } from "@geist-ui/react";

const menu_data = [
  {
    key: 0,
    name: "Characterization",
    publicId: "primerdriver/characterize",
    href: "/characterize",
    color: "error",
  },
  {
    key: 1,
    name: "DNA-based",
    publicId: "primerdriver/dna",
    href: "/dna",
    color: "success",
  },
  {
    key: 2,
    name: "Protein-based",
    publicId: "primerdriver/protein.png",
    href: "/protein",
    color: "cyan",
  },
];

function Menu() {
  return (
    <Grid.Container gap={3} justify="center">
      {menu_data.map((item, i) => (
        <Grid key={i} xs={24} md={6}>
          <Card hoverable>
            <Link id={item.href.slice(1)} to={item.href}>
              <CloudImage
                className="img-fluid"
                cloudName="kdphotography-assets"
                publicId={item.publicId}
                dpr="auto"
                responsive
                responsiveUseBreakpoints
                secure
                width="auto"
                crop="scale"
              />
            </Link>
            <Card.Content>
              <Row justify="center">
                <Link id={item.href.slice(1)} to={item.href}>
                  <Button type={item.color} ghost>
                    {item.name}
                  </Button>
                </Link>
              </Row>
            </Card.Content>
          </Card>
        </Grid>
      ))}
    </Grid.Container>
  );
}

Menu.propTypes = {
  stations: PropTypes.arrayOf(
    PropTypes.shape({
      key: PropTypes.number.isRequired,
      name: PropTypes.string.isRequired,
      publicId: PropTypes.string.isRequired,
      href: PropTypes.string.isRequired,
      color: PropTypes.string.isRequired,
    }),
  ),
};

export default Menu;
