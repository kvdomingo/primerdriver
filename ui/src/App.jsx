import { MDBCol as Col, MDBRow as Row } from "mdbreact";
import "./App.css";
import AppRoutes from "./components/AppRoutes";
import Landing from "./components/Landing";

function App() {
  return (
    <Row className="vh-100 vw-100 mx-0 align-items-center">
      <Col sm={3}>
        <Landing />
      </Col>
      <Col sm={9}>
        <AppRoutes id="app" />
      </Col>
    </Row>
  );
}

export default App;
