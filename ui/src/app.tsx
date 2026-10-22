import { RouterProvider } from "@tanstack/react-router";
import { MDBCol as Col, MDBRow as Row } from "mdbreact";
import "./app.css";
import { router } from "./components/app-routes";
import Landing from "./components/landing";

function App() {
  return (
    <Row className="vh-100 vw-100 mx-0 align-items-center">
      <Col sm="3">
        <Landing />
      </Col>
      <Col sm="9">
        <RouterProvider router={router} />
      </Col>
    </Row>
  );
}

export default App;
