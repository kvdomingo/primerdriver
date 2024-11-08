import { Image } from "cloudinary-react";
import {
  MDBCard as Card,
  MDBCardBody as CardBody,
  MDBCol as Col,
  MDBRow as Row,
} from "mdbreact";
import { useEffect } from "react";
import { Link } from "react-router-dom";
import { usePrimerDriverContext } from "../../contexts/PrimerDriverContext.jsx";

const menu_data = [
  {
    key: 0,
    name: "Characterization",
    publicId: "primerdriver/characterize",
    href: "/characterize",
    color: "danger",
  },
  {
    key: 1,
    name: "DNA-based",
    publicId: "primerdriver/dna",
    href: "/dna",
    color: "primary",
  },
  {
    key: 2,
    name: "Protein-based",
    publicId: "primerdriver/protein.png",
    href: "/protein",
    color: "success",
  },
];

function Menu() {
  const { PDDispatch } = usePrimerDriverContext();

  useEffect(() => {
    PDDispatch({
      type: "updateResults",
      payload: {
        data: [],
        loaded: false,
      },
    });
    PDDispatch({
      type: "updateMode",
      payload: "",
    });
  }, [PDDispatch]);

  return (
    <Row className="row-cols-1 row-cols-md-3 mx-0 h-100 align-items-center">
      {menu_data.map(item => (
        <Col key={item.publicId}>
          <Card>
            <div className="view overlay">
              <Image
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
              <Link to={item.href}>
                <div className="mask rgba-black-slight" />
              </Link>
            </div>
            <CardBody className="text-center">
              <Link
                id={item.href.slice(1)}
                to={item.href}
                className={`btn btn-${item.color} btn-md`}
              >
                {item.name}
              </Link>
            </CardBody>
          </Card>
        </Col>
      ))}
    </Row>
  );
}

export default Menu;
